#!/usr/bin/env python3
"""
Diagnose emergent seasonality of time-varying ERSEM parameters.

For each parameter it calculates:

1. Monthly ensemble-mean value for each year.
2. 12-month climatology in q-space.
3. 12-month climatology converted back to physical ERSEM units.
4. Annual harmonic after accounting for a linear long-term trend:
       q(t) = c + trend*t + a*sin(2*pi*t) + b*cos(2*pi*t)
5. Seasonal amplitude:
       A = sqrt(a^2 + b^2)
6. Partial seasonal R2:
       improvement obtained by adding the annual harmonic
       relative to a trend-only model.
7. Interannual repeatability:
       median correlation between monthly cycles from different years,
       where sufficient overlapping months exist.

Outputs:
    monthly_parameter_climatology.csv
    parameter_seasonality_metrics.csv
    P_parameter_monthly_seasonality.png
    Z_parameter_monthly_seasonality.png
"""

from pathlib import Path
import argparse
import itertools
import re

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


CENTRE = {
    "P1": {"alpha": 4.0, "beta": 0.07, "sum": 1.375, "srs": 0.04},
    "P2": {"alpha": 5.0, "beta": 0.10, "sum": 1.625, "srs": 0.04},
    "P3": {"alpha": 6.0, "beta": 0.12, "sum": 2.0, "srs": 0.045},
    "P4": {"alpha": 3.0, "beta": 0.06, "sum": 1.125, "srs": 0.035},
    "Z4": {"sum": 1.0, "pu": 0.6},
    "Z5": {"sum": 1.25, "pu": 0.5},
    "Z6": {"sum": 1.5, "pu": 0.4},
}

PAT = re.compile(r"^instances_([PZ]\d+)_parameters_(.+)$")

MONTH_NAMES = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
]


def logit(x):
    x = np.clip(x, 1e-12, 1.0 - 1e-12)
    return np.log(x / (1.0 - x))


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def to_q(theta, taxon, param):
    theta0 = CENTRE[taxon][param]

    if param == "pu":
        return logit(theta) - logit(theta0)

    return np.log(np.maximum(theta, 1e-300) / theta0)


def from_q(q, taxon, param):
    theta0 = CENTRE[taxon][param]

    if param == "pu":
        return logistic(logit(theta0) + q)

    return theta0 * np.exp(q)


def read_qmean(ds, analysis_step):
    data = {}

    for var in ds.data_vars:
        m = PAT.match(var)
        if m is None:
            continue

        taxon, param = m.groups()

        if taxon not in CENTRE or param not in CENTRE[taxon]:
            continue

        da = ds[var]

        if "da_step" in da.dims:
            da = da.isel(da_step=analysis_step)

        for dim in list(da.dims):
            if dim not in ("time", "member") and da.sizes[dim] == 1:
                da = da.isel({dim: 0})

        da = da.transpose("time", "member")

        theta = np.asarray(da.values, dtype=float)
        q = to_q(theta, taxon, param)

        data[(taxon, param)] = np.nanmean(q, axis=1)

    return data


def monthly_year_table(times, qmean):
    """
    DEEP: First average all DA cycles within each calendar year-month.

    This prevents a month with more observation profiles from receiving
    disproportionately high weight in the climatology.
    """

    df = pd.DataFrame({
        "time": times,
        "q": qmean,
    })

    df["year"] = df["time"].dt.year
    df["month"] = df["time"].dt.month

    ym = (
        df.groupby(["year", "month"], as_index=False)["q"]
        .mean()
    )

    return ym


def annual_harmonic(times, qmean):
    """
    Compare:
       trend-only model
    against:
       trend + annual sin/cos model.

    Therefore slow drift alone cannot produce a large seasonal score.
    """

    times = pd.DatetimeIndex(times)

    t_days = (
        times - times[0]
    ).total_seconds().to_numpy() / 86400.0

    t_year = t_days / 365.2425

    y = np.asarray(qmean, dtype=float)

    ok = np.isfinite(y) & np.isfinite(t_year)

    y = y[ok]
    t = t_year[ok]

    if len(y) < 10:
        return np.nan, np.nan, np.nan

    Xtrend = np.column_stack([
        np.ones_like(t),
        t,
    ])

    Xseason = np.column_stack([
        np.ones_like(t),
        t,
        np.sin(2.0 * np.pi * t),
        np.cos(2.0 * np.pi * t),
    ])

    b0, *_ = np.linalg.lstsq(Xtrend, y, rcond=None)
    b1, *_ = np.linalg.lstsq(Xseason, y, rcond=None)

    pred0 = Xtrend @ b0
    pred1 = Xseason @ b1

    sse0 = np.sum((y - pred0) ** 2)
    sse1 = np.sum((y - pred1) ** 2)

    if sse0 > 0.0:
        partial_r2 = 1.0 - sse1 / sse0
    else:
        partial_r2 = np.nan

    a = b1[2]
    b = b1[3]

    amplitude = np.sqrt(a*a + b*b)

    # Peak timing of fitted annual component.
    phase = np.arctan2(a, b)

    peak_fraction = (phase / (2.0 * np.pi)) % 1.0
    peak_month = 1.0 + 12.0 * peak_fraction

    return amplitude, partial_r2, peak_month


def repeatability(ym):
    """
    Median pairwise Pearson correlation between individual-year
    monthly cycles.

    Require >= 6 common months for a year pair.
    """

    pivot = ym.pivot(
        index="month",
        columns="year",
        values="q",
    )

    correlations = []

    years = list(pivot.columns)

    for y1, y2 in itertools.combinations(years, 2):
        a = pivot[y1]
        b = pivot[y2]

        ok = a.notna() & b.notna()

        if ok.sum() < 6:
            continue

        aa = a[ok].to_numpy()
        bb = b[ok].to_numpy()

        if np.std(aa) == 0 or np.std(bb) == 0:
            continue

        correlations.append(
            np.corrcoef(aa, bb)[0, 1]
        )

    if not correlations:
        return np.nan, 0

    return float(np.nanmedian(correlations)), len(correlations)


def analyse_parameter(times, qmean, taxon, param):
    ym = monthly_year_table(times, qmean)

    climatology = (
        ym.groupby("month")["q"]
        .agg(["mean", "std", "count"])
        .reindex(range(1, 13))
    )

    climatology.index.name = "month"

    amplitude, partial_r2, peak_month = annual_harmonic(
        times,
        qmean,
    )

    repeat_r, n_pairs = repeatability(ym)

    rows = []

    for month in range(1, 13):
        qmu = climatology.loc[month, "mean"]
        qsd = climatology.loc[month, "std"]
        n = climatology.loc[month, "count"]

        theta = (
            from_q(qmu, taxon, param)
            if np.isfinite(qmu)
            else np.nan
        )

        rows.append({
            "taxon": taxon,
            "parameter": param,
            "month": month,
            "month_name": MONTH_NAMES[month - 1],
            "q_climatology_mean": qmu,
            "q_interannual_sd": qsd,
            "n_year_months": n,
            "physical_parameter": theta,
        })

    metric = {
        "taxon": taxon,
        "parameter": param,
        "annual_harmonic_amplitude_q": amplitude,
        "seasonal_partial_R2": partial_r2,
        "peak_month": peak_month,
        "interannual_repeatability_r": repeat_r,
        "n_year_pairs": n_pairs,
        "q_total_range": np.nanmax(qmean) - np.nanmin(qmean),
    }

    return ym, pd.DataFrame(rows), metric


def plot_group(results, taxa, params, outfile):
    fig, axes = plt.subplots(
        len(taxa),
        len(params),
        figsize=(4.1 * len(params), 2.8 * len(taxa)),
        sharex=True,
        squeeze=False,
    )

    months = np.arange(1, 13)

    for i, taxon in enumerate(taxa):
        for j, param in enumerate(params):

            ax = axes[i, j]
            key = (taxon, param)

            if key not in results:
                ax.axis("off")
                continue

            ym, clim, metric = results[key]

            # Individual-year monthly cycles.
            for year, g in ym.groupby("year"):
                g = g.sort_values("month")
                ax.plot(
                    g["month"],
                    g["q"],
                    linewidth=0.8,
                    alpha=0.25,
                )

            qmu = clim["q_climatology_mean"].to_numpy()
            qsd = clim["q_interannual_sd"].to_numpy()

            ax.fill_between(
                months,
                qmu - qsd,
                qmu + qsd,
                alpha=0.20,
            )

            ax.plot(
                months,
                qmu,
                linewidth=2.0,
                marker="o",
                markersize=3,
            )

            ax.axhline(
                0.0,
                linestyle="--",
                linewidth=0.8,
            )

            r2 = metric["seasonal_partial_R2"]
            rr = metric["interannual_repeatability_r"]
            amp = metric["annual_harmonic_amplitude_q"]

            ax.set_title(
                f"{taxon} {param}\n"
                f"A={amp:.2f}, R²season={r2:.2f}, "
                f"Rrepeat={rr:.2f}"
            )

            ax.set_xlim(1, 12)
            ax.set_xticks(months)
            ax.set_xticklabels(
                MONTH_NAMES,
                rotation=45,
                ha="right",
                fontsize=8,
            )

            ax.grid(alpha=0.2)

            if j == 0:
                ax.set_ylabel("q")

    fig.suptitle(
        "Emergent monthly parameter variability",
        fontsize=14,
    )

    fig.tight_layout()
    fig.savefig(outfile, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        default="out/DAout.nc",
    )

    parser.add_argument(
        "--dir",
        default="vald_seasonality",
    )

    parser.add_argument(
        "--analysis-step",
        type=int,
        default=1,
    )

    args = parser.parse_args()

    outdir = Path(args.dir)
    outdir.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(args.input) as ds:
        times = pd.to_datetime(ds["time"].values)
        qdata = read_qmean(ds, args.analysis_step)

    results = {}
    clim_all = []
    metrics = []

    for (taxon, param), qmean in sorted(qdata.items()):

        ym, clim, metric = analyse_parameter(
            times,
            qmean,
            taxon,
            param,
        )

        results[(taxon, param)] = (
            ym,
            clim,
            metric,
        )

        clim_all.append(clim)
        metrics.append(metric)

    climatology_df = pd.concat(
        clim_all,
        ignore_index=True,
    )

    metrics_df = pd.DataFrame(metrics)

    climatology_df.to_csv(
        outdir / "monthly_parameter_climatology.csv",
        index=False,
    )

    metrics_df.to_csv(
        outdir / "parameter_seasonality_metrics.csv",
        index=False,
    )

    plot_group(
        results,
        taxa=["P1", "P2", "P3", "P4"],
        params=["alpha", "beta", "sum", "srs"],
        outfile=outdir / "P_parameter_monthly_seasonality.png",
    )

    plot_group(
        results,
        taxa=["Z4", "Z5", "Z6"],
        params=["sum", "pu"],
        outfile=outdir / "Z_parameter_monthly_seasonality.png",
    )

    print()
    print("===== DEEP: SEASONALITY METRICS =====")

    show = metrics_df.sort_values(
        ["seasonal_partial_R2",
         "interannual_repeatability_r"],
        ascending=False,
    )

    print(
        show.to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}",
        )
    )

    print()
    print("Done! ;)")
    print(outdir / "monthly_parameter_climatology.csv")
    print(outdir / "parameter_seasonality_metrics.csv")
    print(outdir / "P_parameter_monthly_seasonality.png")
    print(outdir / "Z_parameter_monthly_seasonality.png")


if __name__ == "__main__":
    main()

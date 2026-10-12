#!/usr/bin/env bash
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
set -euo pipefail

# ==============================================================================
# CAMP_L4 generic EAT/PDAF assimilation runner
#
# N_ENS must be supplied explicitly for every experiment.
# Scheduler resources are deliberately not defined here.
# ==============================================================================

: "${N_ENS:?ERROR: set N_ENS explicitly, for example N_ENS=100}"

MPI_PY="${MPI_PY:-1}"
MPI_GOTM="${MPI_GOTM:-${N_ENS}}"

PYTHON_BIN="${PYTHON_BIN:-python}"
MPIEXEC="${MPIEXEC:-mpiexec}"
GOTM_BIN="${GOTM_BIN:-eat-gotm}"

OUTDIR="${OUTDIR:-out}"

# Restart ensemble perturbation
RESTART_LNSIGMA="${RESTART_LNSIGMA:-0.05}"
RESTART_RANDOM_SEED="${RESTART_RANDOM_SEED:-42}"

# Initial parameter-ensemble perturbations
PARAMETER_RANDOM_SEED="${PARAMETER_RANDOM_SEED:-42}"


# Plugin parameter evolution
POS_EPS="${POS_EPS:-1.0e-12}"

export \
    N_ENS \
    RESTART_LNSIGMA \
    RESTART_RANDOM_SEED \
    PARAMETER_RANDOM_SEED \
    POS_EPS

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

for setting in \
    "N_ENS=${N_ENS}" \
    "MPI_PY=${MPI_PY}" \
    "MPI_GOTM=${MPI_GOTM}"
do
    name="${setting%%=*}"
    value="${setting#*=}"

    if ! is_positive_integer "${value}"; then
        echo "ERROR: ${name} must be a positive integer; received '${value}'" >&2
        exit 1
    fi
done

if (( N_ENS > 9999 )); then
    echo "ERROR: N_ENS must not exceed 9999 with four-digit member filenames" >&2
    exit 1
fi

if (( MPI_GOTM != N_ENS )); then
    echo "WARNING: MPI_GOTM=${MPI_GOTM}, but N_ENS=${N_ENS}." >&2
    echo "The tested layout uses one GOTM process per ensemble member." >&2
fi

echo "================================================================================"
echo "CAMP_L4 EAT/PDAF ASSIMILATION"
echo "Working directory:       $(pwd)"
echo "Ensemble members:        ${N_ENS}"
echo "Python MPI tasks:        ${MPI_PY}"
echo "GOTM MPI tasks:          ${MPI_GOTM}"
echo "Restart ln-sigma:        ${RESTART_LNSIGMA}"
echo "Restart random seed:     ${RESTART_RANDOM_SEED}"
echo "Parameter random seed:   ${PARAMETER_RANDOM_SEED}"
echo "Output directory:        ${OUTDIR}"
echo "================================================================================"

for file in \
    gotm.yaml \
    fabm.yaml \
    config/obs.yaml \
    run.py \
    science_state_plugin.py \
    parameter_transform_plugin.py \
    physical_posterior_plugin.py \
    parameter_prior.py \
    parameter_evolution_plugin.py
do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: required file missing or empty: ${file}" >&2
        exit 1
    fi
done

if [[ ! -s "../../spin/restart.nc" ]]; then
    echo "ERROR: missing or empty spin restart: ../../spin/restart.nc" >&2
    exit 1
fi

for command in \
    "${PYTHON_BIN}" \
    "${MPIEXEC}" \
    "${GOTM_BIN}"
do
    if ! command -v "${command}" >/dev/null 2>&1; then
        echo "ERROR: executable not found in PATH: ${command}" >&2
        exit 1
    fi
done

"${PYTHON_BIN}" - <<'PY'
from pathlib import Path
import yaml

gotm_path = Path("gotm.yaml")
fabm_path = Path("fabm.yaml")
obs_path = Path("config/obs.yaml")

gotm = yaml.safe_load(gotm_path.read_text(encoding="utf-8")) or {}
fabm = yaml.safe_load(fabm_path.read_text(encoding="utf-8")) or {}
obs = yaml.safe_load(obs_path.read_text(encoding="utf-8"))

required_gotm = {
    "time.start": gotm.get("time", {}).get("start"),
    "time.stop": gotm.get("time", {}).get("stop"),
    "fabm.yaml_file": gotm.get("fabm", {}).get("yaml_file"),
}

missing = [
    name
    for name, value in required_gotm.items()
    if value in (None, "")
]

if missing:
    raise SystemExit(
        "ERROR: missing required GOTM settings: "
        + ", ".join(missing)
    )

if not fabm.get("instances"):
    raise SystemExit(
        "ERROR: fabm.yaml does not contain an instances mapping"
    )

if not isinstance(obs, dict):
    raise SystemExit(
        "ERROR: config/obs.yaml must contain a top-level mapping"
    )

streams = obs.get("streams")

if not isinstance(streams, list) or not streams:
    raise SystemExit(
        "ERROR: config/obs.yaml must contain a non-empty streams list"
    )

for index, stream in enumerate(streams, start=1):
    if not isinstance(stream, dict):
        raise SystemExit(
            f"ERROR: observation stream {index} is not a mapping"
        )

    missing_stream_keys = [
        key for key in ("var", "file")
        if not stream.get(key)
    ]

    if missing_stream_keys:
        raise SystemExit(
            f"ERROR: observation stream {index} is missing: "
            + ", ".join(missing_stream_keys)
        )

    observation_file = Path(stream["file"])

    if not observation_file.is_file():
        raise SystemExit(
            f"ERROR: observation file does not exist: "
            f"{observation_file}"
        )

    if observation_file.stat().st_size == 0:
        raise SystemExit(
            f"ERROR: observation file is empty: "
            f"{observation_file}"
        )

restart = gotm.get("restart", {})

if restart.get("load") is not True:
    raise SystemExit(
        "ERROR: assimilation gotm.yaml must use restart.load: true"
    )

print("Configuration checks passed")
print("  Start:              ", gotm["time"]["start"])
print("  Stop:               ", gotm["time"]["stop"])
print("  FABM:               ", gotm["fabm"]["yaml_file"])
print("  Observation streams:", len(streams))
PY

# DEEP: Do not overwrite a previous assimilation accidentally.
for pattern in \
    "restart_*.nc" \
    "gotm_*.yaml" \
    "fabm_*.yaml" \
    "gotm_*.stdout" \
    "gotm_*.stderr"
do
    if compgen -G "${pattern}" >/dev/null; then
        echo "ERROR: existing runtime files match '${pattern}'." >&2
        echo "Use a fresh run directory or archive the previous attempt first." >&2
        exit 1
    fi
done

if [[ -s "${OUTDIR}/DAout.nc" ]]; then
    echo "ERROR: existing assimilation output found: ${OUTDIR}/DAout.nc" >&2
    exit 1
fi

mkdir -p "${OUTDIR}" run_logs

cp -f ../../spin/restart.nc restart.nc

echo
echo "[1/3] Creating ${N_ENS} member restart files..."

"${PYTHON_BIN}" - <<'PY'
import os
import shutil

import eatpy
import numpy as np

n_ens = int(os.environ["N_ENS"])
sigma = float(os.environ["RESTART_LNSIGMA"])
seed = int(os.environ["RESTART_RANDOM_SEED"])

if sigma <= 0.0:
    for member in range(1, n_ens + 1):
        shutil.copyfile(
            "restart.nc",
            f"restart_{member:04d}.nc",
        )

    print(
        f"Copied restart.nc to {n_ens} ensemble members "
        "without perturbation"
    )

else:
    rng = np.random.default_rng(seed)

    with eatpy.models.gotm.RestartEnsemble(
        "restart.nc",
        n_ens,
    ) as ensemble:

        for name, values in ensemble.template.items():
            if name in {"h", "z", "zi", "lon", "lat"}:
                continue

            factors = rng.lognormal(
                mean=0.0,
                sigma=sigma,
                size=(n_ens,) + values.shape,
            )

            # Preserve the sign of originally negative restart values.
            # Retain the ARGO draws and floor for non-negative values.
            original = np.asarray(values[:])
            perturbed = original * factors
            ensemble[name] = np.where(
                original < 0.0,
                perturbed,
                np.maximum(perturbed, 1.0e-12),
            )

    print(
        f"Created {n_ens} perturbed member restarts "
        f"using sigma={sigma} and seed={seed}"
    )
PY

echo
echo "[2/3] Creating member-specific GOTM and FABM configurations..."

"${PYTHON_BIN}" - <<'PYCODE'
import os

import eatpy

from parameter_prior import (
    load_centres,
    prior_range,
    draw_initial_prior,
)


n_ens = int(
    os.environ["N_ENS"]
)

seed = int(
    os.environ["PARAMETER_RANDOM_SEED"]
)

import numpy as np

rng = np.random.default_rng(
    seed
)


phytoplankton = (
    "P1",
    "P2",
    "P3",
    "P4",
)

zooplankton = (
    "Z4",
    "Z5",
    "Z6",
)


parameter_specs = (
    [
        (
            f"instances/{group}/parameters/{parameter}",
            group,
            parameter,
        )
        for group in phytoplankton
        for parameter in (
            "alpha",
            "beta",
            "sum",
            "srs",
        )
    ]
    +
    [
        (
            f"instances/{group}/parameters/{parameter}",
            group,
            parameter,
        )
        for group in zooplankton
        for parameter in (
            "sum",
            "pu",
        )
    ]
)


if len(parameter_specs) != 22:
    raise RuntimeError(
        "Expected 22 CAMP candidate parameters; "
        f"found {len(parameter_specs)}"
    )


centres = load_centres(
    "fabm.yaml"
)


def repeat(value):
    return [value] * n_ens


print(
    "Active CAMP candidate parameters:",
    len(parameter_specs),
)

print(
    "Initial parameter prior: "
    "uniform +/-30% around standard ERSEM values"
)

print(
    "IMPORTANT: +/-30% applies to INITIAL PRIOR ONLY; "
    "it is not a persistent DA bound."
)


with (
    eatpy.models.gotm.YAMLEnsemble(
        "gotm.yaml",
        n_ens,
    ) as gotm_ensemble,
    eatpy.models.gotm.YAMLEnsemble(
        "fabm.yaml",
        n_ens,
    ) as fabm_ensemble,
):

    gotm_ensemble[
        "fabm/yaml_file"
    ] = fabm_ensemble.file_paths

    gotm_ensemble[
        "restart/load"
    ] = repeat(True)

    for (
        path,
        group,
        parameter,
    ) in parameter_specs:

        centre = centres[
            (group, parameter)
        ]

        lower, upper = prior_range(
            centre,
            parameter,
        )

        samples = draw_initial_prior(
            rng,
            centre,
            parameter,
            n_ens,
        )

        # Sanity checks only.
        # There is deliberately NO clipping.
        if not np.all(
            np.isfinite(samples)
        ):
            raise RuntimeError(
                f"Non-finite initial prior for {path}"
            )

        if parameter == "pu":

            if not np.all(
                (samples > 0.0)
                & (samples < 1.0)
            ):
                raise RuntimeError(
                    f"Invalid pu initial prior for {path}"
                )

        elif not np.all(
            samples > 0.0
        ):
            raise RuntimeError(
                f"Non-positive initial prior for {path}"
            )

        fabm_ensemble[
            path
        ] = samples

        print(
            f"  {path}: "
            f"centre={centre:.8g}, "
            f"initial_prior=U("
            f"{lower:.8g},{upper:.8g}), "
            f"sample_mean={np.mean(samples):.8g}, "
            f"sample_sd={np.std(samples):.8g}"
        )


print(
    f"Created {n_ens} member configurations "
    f"using parameter seed {seed}"
)
PYCODE

cat > run_configuration.txt <<CONFIG
N_ENS=${N_ENS}
MPI_PY=${MPI_PY}
MPI_GOTM=${MPI_GOTM}
RESTART_LNSIGMA=${RESTART_LNSIGMA}
RESTART_RANDOM_SEED=${RESTART_RANDOM_SEED}
PARAMETER_RANDOM_SEED=${PARAMETER_RANDOM_SEED}
POS_EPS=${POS_EPS}
PARAM_TAU_DAYS=${PARAM_TAU_DAYS}
PARAM_PRIOR_FRACTION=${PARAM_PRIOR_FRACTION}
PARAM_EVOLUTION_SEED=42001
LOG_OBS_SD_MIN=${LOG_OBS_SD_MIN}
INITIAL_PRIOR_FRACTION=0.30
RESTART_NEGATIVE_POLICY=preserve_sign
CONFIG

echo
echo "[3/3] Starting EAT/PDAF assimilation..."

"${MPIEXEC}" \
    -n "${MPI_PY}" \
    "${PYTHON_BIN}" -u run.py \
    : \
    -n "${MPI_GOTM}" \
    "${GOTM_BIN}" \
    --separate_gotm_yaml \
    --separate_restart_file

echo
echo "[post] Collecting generated runtime files..."

shopt -s nullglob

runtime_files=(
    result*.nc
    *_surface_daily_mean*.nc
    *_bottom_daily_mean*.nc
    restart_*.nc
    fabm_*.yaml
    gotm_*.yaml
)

for file in "${runtime_files[@]}"; do
    mv -f "${file}" "${OUTDIR}/"
done

log_files=(
    gotm_*.stdout
    gotm_*.stderr
)

for file in "${log_files[@]}"; do
    mv -f "${file}" run_logs/
done

shopt -u nullglob

if [[ ! -s "${OUTDIR}/DAout.nc" ]]; then
    echo "ERROR: expected output was not produced: ${OUTDIR}/DAout.nc" >&2
    exit 2
fi

echo
echo "Assimilation completed successfully"
echo "DA output: $(readlink -f "${OUTDIR}/DAout.nc")"
echo "Runtime configuration: $(readlink -f run_configuration.txt)"

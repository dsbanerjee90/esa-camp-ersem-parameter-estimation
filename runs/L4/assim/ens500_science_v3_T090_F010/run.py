#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

import inspect
import logging
import os
import yaml

import eatpy

from science_state_plugin import ScienceStatePlugin
from physical_posterior_plugin import PhysicalPosteriorWriter
from parameter_transform_plugin import ParameterTransform
from parameter_evolution_plugin import ParameterEvolution


logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

OUTFILE = "out/DAout.nc"
OBS_CFG = "config/obs.yaml"


# ============================================================
# CAMP candidate parameter set
# ============================================================

P = ("P1", "P2", "P3", "P4")
Z = ("Z4", "Z5", "Z6")

PARAMS_SLASH = (
    [
        f"instances/{p}/parameters/{k}"
        for p in P
        for k in ("alpha", "beta", "sum", "srs")
    ]
    +
    [
        f"instances/{z}/parameters/{k}"
        for z in Z
        for k in ("sum", "pu")
    ]
)

if len(PARAMS_SLASH) != 22:
    raise RuntimeError(
        f"Expected 22 CAMP candidate parameters; "
        f"found {len(PARAMS_SLASH)}"
    )


def to_underscore(paths):
    return [
        path.replace("/", "_")
        for path in paths
    ]


PARAMS_UNDERSCORE = to_underscore(
    PARAMS_SLASH
)


def main():

    # ========================================================
    # GOTM / ERSEM augmented state
    # ========================================================

    gotm_kwargs = dict(
        diagnostics_in_state=[
            "total_chlorophyll",
            "O2_o",
            "N3_n",
            "N1_p",
            "N5_s",
            "temp",
            "salt",
        ]
    )

    sig = inspect.signature(
        eatpy.models.GOTM
    )

    if "fabm_parameters_in_state" in sig.parameters:

        gotm_kwargs[
            "fabm_parameters_in_state"
        ] = PARAMS_SLASH

        LOG.info(
            "GOTM(): using fabm_parameters_in_state "
            "with %d parameters.",
            len(PARAMS_SLASH),
        )

    elif "parameters_in_state" in sig.parameters:

        gotm_kwargs[
            "parameters_in_state"
        ] = PARAMS_UNDERSCORE

        LOG.info(
            "GOTM(): using parameters_in_state "
            "with %d parameters.",
            len(PARAMS_UNDERSCORE),
        )

    else:

        raise RuntimeError(
            "Installed EAT/GOTM exposes no supported "
            "FABM parameter-in-state interface."
        )

    experiment = eatpy.models.GOTM(
        **gotm_kwargs
    )

    # ========================================================
    # Immediately before PDAF:
    #
    #   chlorophyll is in log10 space
    #   parameters are in q-space
    #
    # AFTER analysis executes in REVERSE registration order:
    #
    #   ParameterTransform
    #       q -> physical theta
    #
    #   PhysicalPosteriorWriter
    #       save physical observation-constrained posterior
    #
    #   ScienceStatePlugin
    #       chlorophyll -> physical concentration
    #       physical-state safeguards
    #
    #   NetCDF
    #       save final physical analysis state
    # ========================================================

    experiment.add_plugin(
        eatpy.plugins.output.NetCDF(
            OUTFILE
        )
    )

    experiment.add_plugin(
        ScienceStatePlugin()
    )

    experiment.add_plugin(
        PhysicalPosteriorWriter()
    )

    experiment.add_plugin(
        ParameterTransform(
            "fabm.yaml"
        )
    )

    # ========================================================
    # Generic parameter evolution in transformed q-space.
    #
    # Registered AFTER ParameterTransform deliberately:
    #
    # BEFORE analysis:
    #   ParameterTransform: physical theta -> q
    #   ParameterEvolution: records current assimilation time
    #
    # AFTER analysis executes in reverse:
    #   ParameterEvolution: evolves analysed q for next forecast
    #   ParameterTransform: q -> physical theta
    # ========================================================

    tau_days = float(
        os.environ["PARAM_TAU_DAYS"]
    )

    prior_fraction = float(
        os.environ["PARAM_PRIOR_FRACTION"]
    )

    experiment.add_plugin(
        ParameterEvolution(
            obs_yaml="config/obs.yaml",
            fabm_yaml="fabm.yaml",
            tau_days=tau_days,
            prior_fraction=prior_fraction,
            seed=42001,
        )
    )

    LOG.info(
        "Parameter evolution enabled: "
        "tau_days=%g prior_fraction=%g",
        tau_days,
        prior_fraction,
    )

    # ========================================================
    # PDAF
    # ========================================================

    filt = eatpy.PDAF(
        eatpy.pdaf.FilterType.ESTKF
    )

    # ========================================================
    # Observations
    # ========================================================

    with open(
        OBS_CFG,
        "r",
        encoding="utf-8",
    ) as f:

        cfg = yaml.safe_load(f) or {}

    streams = cfg.get(
        "streams",
        [],
    )

    if not streams:
        raise RuntimeError(
            "No observation streams found in "
            f"{OBS_CFG}"
        )

    for stream in streams:

        experiment.add_observations(
            stream["var"],
            stream["file"],
        )

    experiment.run(
        filt
    )


if __name__ == "__main__":
    main()

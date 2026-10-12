#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

from typing import List
import logging
import os
import re

import numpy as np
import eatpy.shared


LOG = logging.getLogger("science_state_plugin")


# ============================================================
# State variables retained in the DA state
# ============================================================

BGC_VARS: List[str] = [
    "total_chlorophyll",
    "O2_o",
    "N3_n",
    "N1_p",
    "N5_s",
    "temp",
    "salt",
]

LOG_VARS: List[str] = [
    "total_chlorophyll",
]

SMALL = float(
    os.getenv("POS_EPS", "1.0e-12")
)

# DEEP:  Minimum chlorophyll observation uncertainty after
# conversion to log10 space.
LOG_OBS_SD_MIN = float(
    os.getenv("LOG_OBS_SD_MIN", "0.05")
)


# DEEP: Exactly the CAMP parameter names that are allowed to remain
# in the filter state.
PARAM_SLASH = re.compile(
    r"^instances/(P[1-4]|Z[4-6])/parameters/"
    r"(alpha|beta|sum|srs|pu)$"
)

PARAM_UNDERSCORE = re.compile(
    r"^instances_(P[1-4]|Z[4-6])_parameters_"
    r"(alpha|beta|sum|srs|pu)$"
)


def _start_stop(info):
    start = int(info["start"])

    stop = int(
        info.get(
            "stop",
            start + int(info["length"]),
        )
    )

    return start, stop


class ScienceStatePlugin(eatpy.shared.Plugin):
    """
    DEEP: State-only plugin for the CAMP science-V2 experiment.

    Responsibilities
    ----------------
    1. Retain the selected observed/model BGC state variables.
    2. Retain the 22 CAMP ERSEM parameters in the filter state.
    3. Transform chlorophyll state and observations to log10 space
       immediately before PDAF.
    4. Back-transform chlorophyll immediately after PDAF.
    5. Prevent physically impossible negative analysed BGC states.

    Unlike my previous versions here the followings are NOT performed
    ------------------------------
    * parameter clipping
    * parameter bounds
    * parameter jitter
    * random walk
    * anti-collapse redistribution
    * parameter shrinkage
    * parameter perturbation

    """

    def initialize(
        self,
        variables,
        ensemble_size,
    ):

        self.vars = variables

        param_names = [
            name
            for name in variables
            if (
                PARAM_SLASH.match(name)
                or PARAM_UNDERSCORE.match(name)
            )
        ]

        if len(param_names) != 22:
            raise RuntimeError(
                "Uh oh! ScienceStatePlugin expected exactly "
                f"22 CAMP parameters; found {len(param_names)}"
            )

        # Preserve exactly the BGC variables required by this experiment
        # plus the 22 parameter state elements.
        keep = (
            set(BGC_VARS)
            | set(param_names)
        )

        for name in list(variables):

            if name not in keep:
                del variables[name]

        self.log_infos = [
            variables[name]
            for name in LOG_VARS
            if name in variables
        ]

        self.param_infos = {
            name: variables[name]
            for name in param_names
        }

        LOG.info(
            "ScienceStatePlugin initialized: "
            "%d CAMP parameters retained; "
            "%d BGC variables requested.",
            len(self.param_infos),
            len(BGC_VARS),
        )

    def before_analysis(
        self,
        time,
        state,
        iobs,
        obs,
        obs_sds,
        filter,
    ):

        # ------------------------------------------------------------
        # DEEP: Chlorophyll forecast state:
        # linear concentration -> log10 concentration (Before assimilation)
        # ------------------------------------------------------------

        for info in self.log_infos:

            linear_state = np.asarray(
                info["data"],
                dtype=float,
            ).copy()

            finite = np.isfinite(
                linear_state
            )

            if not np.all(finite):

                bad = np.argwhere(
                    ~finite
                )

                LOG.error(
                    "Uh oh! Non-finite chlorophyll forecast state: "
                    "count=%d, shape=%s, first_indices=%s",
                    int(bad.shape[0]),
                    tuple(linear_state.shape),
                    bad[:20].tolist(),
                )

                raise RuntimeError(
                    "Uh oh! Non-finite chlorophyll forecast state "
                    "before log10 transformation."
                )

            info["data"][...] = np.log10(
                np.maximum(
                    linear_state,
                    SMALL,
                )
            )

            if not np.all(
                np.isfinite(info["data"])
            ):
                raise RuntimeError(
                    "Uh oh! Non-finite chlorophyll forecast state "
                    "after log10 transformation."
                )

            # --------------------------------------------------------
            # Identify chlorophyll observations mapped to this state
            # --------------------------------------------------------

            start, stop = _start_stop(
                info
            )

            affected = (
                (iobs >= start)
                & (iobs < stop)
            )

            if not np.any(
                affected
            ):
                continue

            linear_obs = np.asarray(
                obs[affected],
                dtype=float,
            ).copy()

            linear_sd = np.asarray(
                obs_sds[affected],
                dtype=float,
            ).copy()

            if (
                not np.all(
                    np.isfinite(linear_obs)
                )
                or not np.all(
                    np.isfinite(linear_sd)
                )
                or np.any(
                    linear_obs <= 0.0
                )
                or np.any(
                    linear_sd <= 0.0
                )
            ):
                raise RuntimeError(
                    "Uh oh! Invalid chlorophyll observation "
                    "or uncertainty before log10 transformation."
                )

            # Delta-method uncertainty transformation:
            #
            # sigma_log10 =
            # sigma_linear / (value * ln(10))
            #
            obs[affected] = np.log10(
                linear_obs
            )

            obs_sds[affected] = np.maximum(
                linear_sd
                / (
                    linear_obs
                    * np.log(10.0)
                ),
                LOG_OBS_SD_MIN,
            )

            if (
                not np.all(
                    np.isfinite(
                        obs[affected]
                    )
                )
                or not np.all(
                    np.isfinite(
                        obs_sds[affected]
                    )
                )
                or np.any(
                    obs_sds[affected]
                    <= 0.0
                )
            ):
                raise RuntimeError(
                    "Uh oh! Invalid chlorophyll observation "
                    "or uncertainty after log10 transformation."
                )


    def after_analysis(
        self,
        state,
    ):

        # ------------------------------------------------------------
        # Chlorophyll:
        # log10 analysis -> physical concentration
        # ------------------------------------------------------------

        for info in self.log_infos:

            info["data"][...] = (
                10.0
                ** info["data"]
            )

        # ------------------------------------------------------------
        # Physical-state safeguards only.
        #
        # These are MODEL STATE VARIABLES, not parameters.
        # ------------------------------------------------------------

        constrained_states = (
            "total_chlorophyll",
            "O2_o",
            "N3_n",
            "N1_p",
            "N5_s",
            "salt",
        )

        for name in constrained_states:

            if name not in self.vars:
                continue

            data = self.vars[name]["data"]

            if not np.all(
                np.isfinite(data)
            ):

                bad = np.argwhere(
                    ~np.isfinite(data)
                )

                raise RuntimeError(
                        f"Oh no! :( Non-finite analysed state for {name}: "
                    f"count={bad.shape[0]}, "
                    f"first_indices={bad[:20].tolist()}"
                )

            n_clipped = int(
                np.count_nonzero(
                    data < SMALL
                )
            )

            if n_clipped:

                LOG.warning(
                    "DEEP: Post-analysis MODEL-STATE constraint: "
                    "%s clipped %d values below %g",
                    name,
                    n_clipped,
                    SMALL,
                )

                np.maximum(
                    data,
                    SMALL,
                    out=data,
                )

        # Temperature may legitimately be below zero.
        if "temp" in self.vars:

            temperature = (
                self.vars["temp"]["data"]
            )

            if not np.all(
                np.isfinite(
                    temperature
                )
            ):

                bad = np.argwhere(
                    ~np.isfinite(
                        temperature
                    )
                )

                raise RuntimeError(
                    "Oh no! Non-finite analysed temperature: "
                    f"count={bad.shape[0]}, "
                    f"first_indices={bad[:20].tolist()}"
                )


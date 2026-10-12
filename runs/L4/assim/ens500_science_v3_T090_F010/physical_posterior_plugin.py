#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

from pathlib import Path
import logging

import numpy as np
import eatpy


LOG = logging.getLogger(__name__)

OUTDIR = Path(
    "out/physical_parameter_posterior"
)


class PhysicalPosteriorWriter(eatpy.shared.Plugin):
    """
    DEEP: Read-only diagnostic writer for CAMP science-V2.

    Intended callback ordering:

        PDAF
          |
          v
        ParameterTransform.after_analysis
          q -> physical theta
          |
          v
        PhysicalPosteriorWriter.after_analysis
          save physical theta
          |
          v
        ScienceStatePlugin.after_analysis
          |
          v
        NetCDF.after_analysis

    This plugin must never modify model/filter state.
    """

    def initialize(
        self,
        variables,
        ensemble_size,
    ):

        self.vars = variables
        self.step = 0
        self.analysis_time = None

        self.param_names = sorted(
            [
                name
                for name in variables
                if (
                    name.startswith("instances/")
                    and "/parameters/" in name
                )
                or (
                    name.startswith("instances_")
                    and "_parameters_" in name
                )
            ]
        )

        if len(self.param_names) != 22:

            raise RuntimeError(
                "PhysicalPosteriorWriter expected "
                f"22 parameter variables; found "
                f"{len(self.param_names)}:\n"
                + "\n".join(
                    self.param_names
                )
            )

        OUTDIR.mkdir(
            parents=True,
            exist_ok=True,
        )

        LOG.info(
            "PhysicalPosteriorWriter initialized: "
            "%d parameters, ensemble_size=%d",
            len(self.param_names),
            ensemble_size,
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

        # Record analysis time only.
        # Do not alter state, observations or uncertainties.
        self.analysis_time = time


    def after_analysis(
        self,
        state,
    ):

        arrays = {}

        for name in self.param_names:

            x = np.asarray(
                self.vars[name]["data"],
                dtype=float,
            )

            if not np.all(
                np.isfinite(x)
            ):
                raise RuntimeError(
                    "Non-finite physical posterior "
                    f"for parameter {name}"
                )

            # DEEP: Copy explicitly so the saved data cannot retain a mutable
            # reference to EAT's filter/model state.
            arrays[name] = x.copy()

        fn = OUTDIR / (
            f"physical_posterior_"
            f"{self.step:04d}.npz"
        )

        np.savez(
            fn,
            analysis_time=str(
                self.analysis_time
            ),
            **arrays,
        )

        LOG.info(
            "Saved physical parameter posterior: "
            "step=%d time=%s file=%s",
            self.step,
            self.analysis_time,
            fn,
        )

        self.step += 1

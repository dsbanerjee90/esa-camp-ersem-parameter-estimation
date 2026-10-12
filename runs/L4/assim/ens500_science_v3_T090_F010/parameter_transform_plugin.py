#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

import logging
import re
from pathlib import Path

import numpy as np
import yaml
import eatpy


LOG = logging.getLogger(__name__)

PARAM_SLASH = re.compile(
    r"^instances/(P[1-4]|Z[4-6])/parameters/"
    r"(alpha|beta|sum|srs|pu)$"
)

PARAM_UNDERSCORE = re.compile(
    r"^instances_(P[1-4]|Z[4-6])_parameters_"
    r"(alpha|beta|sum|srs|pu)$"
)


def _parse_parameter(name):
    """
    Return (group, parameter_type) for one of the 22 CAMP parameters.
    """
    m = PARAM_SLASH.match(name)

    if m is None:
        m = PARAM_UNDERSCORE.match(name)

    if m is None:
        return None

    return m.group(1), m.group(2)


def _load_centres(path="fabm.yaml"):
    """
    DEEP: Load the standard/static ERSEM parameter values from the
    experiment's own fabm.yaml.

    These values define q=0.
    """
    fab = yaml.safe_load(
        Path(path).read_text(encoding="utf-8")
    ) or {}

    instances = fab.get("instances", {})

    centres = {}

    for group in ("P1", "P2", "P3", "P4"):

        for ptype in ("alpha", "beta", "sum", "srs"):

            value = float(
                instances[group]["parameters"][ptype]
            )

            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(
                    f"Invalid standard ERSEM value "
                    f"{group}.{ptype}={value}"
                )

            centres[(group, ptype)] = value

    for group in ("Z4", "Z5", "Z6"):

        for ptype in ("sum", "pu"):

            value = float(
                instances[group]["parameters"][ptype]
            )

            if not np.isfinite(value):
                raise ValueError(
                    f"Invalid standard ERSEM value "
                    f"{group}.{ptype}={value}"
                )

            if ptype == "sum" and value <= 0.0:
                raise ValueError(
                    f"{group}.sum must be > 0"
                )

            if ptype == "pu" and not (0.0 < value < 1.0):
                raise ValueError(
                    f"{group}.pu must be between 0 and 1"
                )

            centres[(group, ptype)] = value

    if len(centres) != 22:
        raise RuntimeError(
            f"Uh oh! Expected 22 parameter centres, got {len(centres)}"
        )

    return centres


def _logit(x):
    return np.log(x / (1.0 - x))


def _inv_logit(x):
    # numerically stable logistic
    x = np.asarray(x, dtype=float)

    out = np.empty_like(x)

    positive = x >= 0.0

    out[positive] = (
        1.0 /
        (1.0 + np.exp(-x[positive]))
    )

    e = np.exp(x[~positive])

    out[~positive] = e / (1.0 + e)

    return out


class ParameterTransform(eatpy.shared.Plugin):
    """
    DEEP: Transform the 22 CAMP ERSEM parameters into coordinates suitable
    for unconstrained ensemble DA.

    Positive parameters:
        q = ln(theta / theta0)

    where theta0 is the standard value in fabm.yaml.

    Assimilation efficiency pu:
        q = logit(theta) - logit(theta0)

    Thus q=0 always corresponds exactly to the standard ERSEM value.

    """

    def __init__(self, fabm_yaml="fabm.yaml"):
        super().__init__()
        self.fabm_yaml = fabm_yaml
        self.centres = _load_centres(fabm_yaml)
        self.parameters = {}

    def initialize(self, variables, ensemble_size):

        for name, info in variables.items():

            parsed = _parse_parameter(name)

            if parsed is None:
                continue

            group, ptype = parsed

            self.parameters[name] = {
                "info": info,
                "group": group,
                "ptype": ptype,
                "centre": self.centres[(group, ptype)],
            }

        if len(self.parameters) != 22:
            raise RuntimeError(
                "Uh oh! ParameterTransform expected exactly 22 CAMP "
                f"parameters in filter state, found "
                f"{len(self.parameters)}"
            )

        LOG.info(
            "DEEP: ParameterTransform initialized for %d parameters",
            len(self.parameters),
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

        for name, meta in self.parameters.items():

            x = meta["info"]["data"]
            ptype = meta["ptype"]
            centre = meta["centre"]

            if not np.all(np.isfinite(x)):
                raise RuntimeError(
                    f"Oh no! Non-finite parameter before transform: {name}"
                )

            if ptype == "pu":

                if np.any((x <= 0.0) | (x >= 1.0)):
                    raise RuntimeError(
                        f"{name}: pu must satisfy 0 < pu < 1 "
                        "before transformation"
                    )

                x[...] = (
                    _logit(x) -
                    _logit(centre)
                )

            else:

                if np.any(x <= 0.0):
                    raise RuntimeError(
                        f"{name}: positive ERSEM parameter "
                        "became <= 0 before transformation"
                    )

                x[...] = np.log(
                    x / centre
                )

    def after_analysis(self, state):

        for name, meta in self.parameters.items():

            q = meta["info"]["data"]
            ptype = meta["ptype"]
            centre = meta["centre"]

            if not np.all(np.isfinite(q)):
                raise RuntimeError(
                    f"Non-finite transformed parameter "
                    f"after analysis: {name}"
                )

            if ptype == "pu":

                q0 = _logit(centre)

                q[...] = _inv_logit(
                    q + q0
                )

            else:

                q[...] = (
                    centre * np.exp(q)
                )

            if not np.all(np.isfinite(q)):
                raise RuntimeError(
                    f"Non-finite physical parameter "
                    f"after inverse transform: {name}"
                )

#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

from pathlib import Path

import numpy as np
import yaml


GROUPS_PHYTO = ("P1", "P2", "P3", "P4")
GROUPS_ZOO = ("Z4", "Z5", "Z6")

PRIOR_FRACTION = 0.30


def load_centres(path="fabm.yaml"):

    fab = yaml.safe_load(
        Path(path).read_text(encoding="utf-8")
    ) or {}

    instances = fab["instances"]

    centres = {}

    for group in GROUPS_PHYTO:
        for parameter in (
            "alpha",
            "beta",
            "sum",
            "srs",
        ):
            value = float(
                instances[group]["parameters"][parameter]
            )

            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(
                    f"Invalid centre {group}.{parameter}={value}"
                )

            centres[(group, parameter)] = value

    for group in GROUPS_ZOO:
        for parameter in ("sum", "pu"):

            value = float(
                instances[group]["parameters"][parameter]
            )

            if parameter == "sum":
                valid = (
                    np.isfinite(value)
                    and value > 0.0
                )
            else:
                valid = (
                    np.isfinite(value)
                    and 0.0 < value < 1.0
                )

            if not valid:
                raise ValueError(
                    f"Invalid centre {group}.{parameter}={value}"
                )

            centres[(group, parameter)] = value

    if len(centres) != 22:
        raise RuntimeError(
            f"Expected 22 parameter centres; got {len(centres)}"
        )

    return centres


def prior_range(
    centre,
    parameter,
    fraction=PRIOR_FRACTION,
):

    if not 0.0 < fraction < 1.0:
        raise ValueError(
            "fraction must satisfy 0 < fraction < 1"
        )

    lower = (
        (1.0 - fraction)
        * centre
    )

    upper = (
        (1.0 + fraction)
        * centre
    )

    if parameter == "pu":

        # This is a sanity assertion, not clipping.
        # Our current ERSEM pu defaults ±30% lie strictly
        # inside the physically meaningful interval (0,1).
        if not (
            0.0 < lower < upper < 1.0
        ):
            raise ValueError(
                f"Requested pu prior [{lower}, {upper}] "
                "is not strictly within (0,1)"
            )

    elif not (
        0.0 < lower < upper
    ):
        raise ValueError(
            f"Invalid positive-parameter prior "
            f"[{lower}, {upper}]"
        )

    return lower, upper


def draw_initial_prior(
    rng,
    centre,
    parameter,
    n_ens,
    fraction=PRIOR_FRACTION,
):

    lower, upper = prior_range(
        centre,
        parameter,
        fraction,
    )

    # IMPORTANT:
    # Draw directly from the prescribed prior.
    return rng.uniform(
        lower,
        upper,
        size=n_ens,
    )

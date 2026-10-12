# ESA-CAMP: Ensemble Data Assimilation and Parameter Estimation in GOTM–ERSEM

**Research software and reproducibility package**
**Author:** Deep S. Banerjee
**Institution:** Plymouth Marine Laboratory (PML), United Kingdom
**Funding:** European Space Agency (ESA), CAMP project

## Overview

This repository provides the code, configurations, input files, and selected diagnostics for a one-dimensional (1-D) GOTM–FABM–ERSEM ensemble data-assimilation and ecosystem-parameter-estimation experiment. The workflow combines the Ensemble and Assimilation Tool (EAT) and the Parallel Data Assimilation Framework (PDAF) with parameter-space transformations, stochastic parameter evolution, and diagnostic analysis.

The packaged L4 experiment demonstrates how biogeochemical observations can be used to constrain time-varying ERSEM parameters. It is a documented research configuration, **not** a standalone distribution of GOTM, FABM, ERSEM, EAT, or PDAF, and is not an operational forecasting product.

## Scientific configuration

| Item | Supplied configuration |
|---|---|
| Test case | L4, 1-D water column |
| Water depth | 50 m |
| Vertical levels | 100 |
| Model time step | 900 s |
| Deterministic spin-up | 2003-02-01 00:30 to 2011-02-01 00:30 |
| Assimilation setup | 2011-02-01 00:30 to 2019-02-01 00:30 |
| Ensemble size | 500 |
| Estimated parameters | 22 ERSEM parameters |
| Assimilation algorithm | Ensemble Square Root Kalman Filter (ESTKF) |
| Initial parameter prior | Uniform range of ±30% around local FABM reference values |
| Parameter memory timescale | 90 days (`T090`) |
| Parameter prior-fraction setting | 0.10 (`F010`) |
| Parameter evolution random seed | 42001 |
| Restart and initial-parameter seeds | 42 |
| Restart perturbation log-standard deviation | 0.05 |

The repository includes a completed spin-up state and an assimilation **setup**; it does not claim to contain a verified complete set of 500-member posterior results.

### Data assimilation and parameter transformations

The EAT/PDAF-based workflow uses ESTKF to update the model state and estimated parameters. The configuration assimilates chlorophyll on a base-10 logarithmic scale. Parameter-space transformations are applied to support physically meaningful constraints, including log-ratio coordinates for positive parameters and a centred-logit representation for the bounded `pu` parameter. Refer to `parameter_prior.py`, `parameter_transform_plugin.py`, and `run.py` for the executable implementation.

### Stochastic parameter evolution

The supplied parameter-evolution setup uses a 90-day timescale and a prior-fraction setting of 0.10. Its conceptual autoregressive form in centred parameter coordinates is:

\[
q(t+\Delta t)=\rho q(t)+\eta,\qquad \rho=\exp(-\Delta t/\tau),
\]

\[
\eta=\sigma_{\infty}\sqrt{1-\rho^2}\,\varepsilon,\qquad \varepsilon\sim\mathcal N(0,1),
\]

where \(\tau=90\) days is the memory timescale and \(\eta\) is the stochastic innovation. The parameterisation and scaling used in an actual run are defined by `parameter_evolution_plugin.py`; the equations here are explanatory rather than a substitute for the source implementation. The `T090_F010` string in the experiment folder is an internal configuration identifier, not a public release version.

The physical-posterior files produced by this workflow may represent **post-evolution** parameters at intermediate schedule times rather than untouched post-analysis values. Users should consult `physical_posterior_plugin.py` before interpreting these products.

## Repository contents

```text
.
├── README.md
├── README_TECHNICAL.md
├── COPY_PROVENANCE.json
├── MANIFEST_SHA256_ORIGINAL.txt
├── MANIFEST_SHA256_RELEASE.txt
├── prepare_spin_rerun.py
├── validate_package.py
├── environment/
│   ├── environment.yml
│   ├── pip-freeze.txt
│   └── runtime_identity.json
└── runs/
    └── L4/
        ├── spin/                         # Spin configuration, forcing, relaxation and restart
        └── assim/
            └── ens500_science_v3_T090_F010/
                ├── run.py
                ├── preflight_l4.py
                ├── parameter_prior.py
                ├── parameter_transform_plugin.py
                ├── parameter_evolution_plugin.py
                ├── physical_posterior_plugin.py
                ├── science_state_plugin.py
                ├── validation_param_seasonality.py
                ├── gotm.yaml
                ├── fabm.yaml
                ├── config/obs.yaml
                ├── obs/
                ├── forcings/
                ├── relax/
                └── vald_seasonality/
```

The inventory above is illustrative. Some other configuration and scheduler files are also supplied; the release manifest is the authoritative list of archived files.

## Observations and model forcing

The assimilation configuration references L4 observations of temperature, salinity, chlorophyll, nitrate, phosphate, silicate and oxygen, together with meteorological forcing and WOA-based relaxation inputs. These external input datasets are not distributed in this GitHub repository. They must be obtained or prepared separately, subject to the relevant data access and licensing conditions.

## Environment and execution

This workflow was prepared in a Linux HPC/Slurm environment with compatible EAT/GOTM/FABM/ERSEM and PDAF installations. The `environment/` folder contains dependency and executable-identity records. It is not a guarantee that the complete model stack can be reproduced solely from the Conda environment specification.

From the repository root, inspect the configuration, paths and Slurm resource request before running any experiment:

```bash
cd runs/L4/assim/ens500_science_v3_T090_F010
python preflight_l4.py
```

Only after preflight passes and resource requirements are approved, submit on a compatible Slurm system:

```bash
sbatch submit.slurm
```

The assimilation setup generates ensemble-specific configuration/restart files. **Do not submit this 500-member experiment without checking allocation, partition, walltime, process layout and storage requirements.**

A completed spin-up state is included. For optional spin-up regeneration without changing the supplied spin results, see `README_TECHNICAL.md` and `prepare_spin_rerun.py`.

## Parameter-seasonality diagnostics

`validation_param_seasonality.py` is included alongside the following existing products under `vald_seasonality/`:

- `monthly_parameter_climatology.csv`
- `parameter_seasonality_metrics.csv`
- `P_parameter_monthly_seasonality.png`
- `Z_parameter_monthly_seasonality.png`

These are provided as example diagnostics. Interpretation should follow the definitions and selection rules in the validation script, not inferred from the figure filenames alone.

## Reproducibility and provenance

The package retains `COPY_PROVENANCE.json` and a historical checksum manifest, `MANIFEST_SHA256_ORIGINAL.txt`, describing an earlier package state. The edited source files in this release are authoritative; the historical manifest is retained as a provenance record, **not** as the integrity check for the release. `MANIFEST_SHA256_RELEASE.txt` should be regenerated after all publication documentation and files have been finalised, and then verified against the exact staged release contents.

Random seeds and environment identifiers are recorded to aid reproducibility. Identical seeds do not guarantee bitwise-identical ensemble members when ensemble size, random-draw shapes, numerical libraries or parallel layout change. Local absolute paths in runtime records may require adjustment and privacy review before redistribution.

## Upstream frameworks and acknowledgements

This work uses the **General Ocean Turbulence Model (GOTM)**, the **Framework for Aquatic Biogeochemical Models (FABM)**, the **European Regional Seas Ecosystem Model (ERSEM)**, the **Ensemble and Assimilation Tool (EAT)** and the **Parallel Data Assimilation Framework (PDAF)**. We acknowledge their original developers and maintainers, particularly Lars Nerger and collaborators for PDAF. The workflow's original research scripts should not be confused with ownership of these upstream frameworks.

PDAF reference:

> Nerger, L., & Hiller, W. (2013). Software for ensemble-based data assimilation systems—Implementation strategies and scalability. *Computers & Geosciences*, **55**, 110–118.

PDAF software: https://doi.org/10.5281/zenodo.7861829

Users should additionally cite the versions and primary scholarly references of EAT, GOTM, FABM and ERSEM actually used in their work. Software redistribution and derivative-work obligations remain governed by each component's applicable licence.

## Funding and attribution

Developed by **Deep S. Banerjee**, Plymouth Marine Laboratory (PML), United Kingdom, with funding from the **European Space Agency (ESA)** through the CAMP project. The formal ESA project identifier and any required institutional acknowledgements should be added once verified.

## Licence, citation and archive

The licence for the original repository components, copyright ownership and third-party redistribution permissions are being verified before the public release. A repository-level `LICENSE`, machine-readable `CITATION.cff`, and a version-specific Zenodo DOI will be added at publication. **Until then, no open-source licence is implied by this README.**

The eventual manuscript Code Availability statement should cite the precise archived release DOI. For technical setup details, refer to `README_TECHNICAL.md`.

## Contact

**Deep S. Banerjee**
Plymouth Marine Laboratory, United Kingdom
GitHub: https://github.com/dsbanerjee90

## ESA CAMP project

This software was developed at Plymouth Marine Laboratory (PML) within the European Space Agency (ESA) CAMP project.

**Official project webpage:** https://climate.esa.int/en/observations-and-modelling/camp/

## Availability of experimental inputs

This GitHub repository distributes the parameter-estimation software,
model configurations, execution scripts, initial restart states,
and selected parameter-seasonality diagnostics.

The original L4 observation streams, meteorological forcing files,
climatological relaxation inputs, and large spin-up diagnostic
outputs are not included in the public GitHub release.

Consequently, the supplied configuration is not independently executable
without obtaining or preparing the required external input datasets.

The original experimental inputs are retained separately in the
ESA-CAMP computational workspace. Users wishing to reproduce the
experiment must provide compatible input data and the required
GOTM–FABM–ERSEM, EAT and PDAF software dependencies.

The included scripts and configuration files document the scientific
implementation, but successful reproduction may require adaptation
to the target computational environment.

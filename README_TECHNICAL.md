# L4 GOTM–ERSEM spin and Science V3: 500 members

Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
Funded by the European Space Agency (ESA).

This package contains the completed L4 spin directory (including
its restart and result NetCDFs), forcing and relaxation inputs, plus a fresh
500-member assimilation setup. No symlinks or hardlinks are created. The source
experiments are not modified and their assimilation runtime results are not copied.

## Science settings

- Spin: 2003-02-01 00:30 to 2011-02-01 00:30; one deterministic GOTM process.
- Assimilation: 2011-02-01 00:30 to 2019-02-01 00:30; 500 members.
- L4: 50 m, 100 levels, dt=900 seconds.
- 22 ERSEM parameters, initial uniform prior +/-30% around local FABM values.
- ESTKF; log10 chlorophyll; log-ratio positive parameters and centred-logit pu.
- Parameter evolution: tau=90 days, prior_fraction=0.10, seed=42001.
- Restart ln-sigma=0.05; restart and initial parameter seeds=42.
- Same V3 plugin order and sign-preserving restart perturbation as the new L4 run.
- physical_parameter_posterior files contain post-evolution physical parameters
  except at the final schedule time, rather than untouched analysis posteriors.
- Same seeds with a different ensemble size do NOT guarantee identical first
  500 members of the 1500-member experiment: draw shapes and normalization differ.

## Check and submit

From the package root, with the eat environment active:

    cd runs/L4/assim/ens500_science_v3_T090_F010
    python preflight_l4.py
    sbatch submit.slurm

The job runs preflight and generates its own 500 restart/configuration files.

## Spin-up

The completed spin restart is already copied. You can start assimilation directly.
To run spin again without altering the supplied results from the package root:

    python prepare_spin_rerun.py
    cd runs/L4/spin_rerun
    sbatch submit.slurm

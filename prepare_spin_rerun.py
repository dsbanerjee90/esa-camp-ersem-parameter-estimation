# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
from pathlib import Path
import shutil
root = Path(__file__).resolve().parent
source = root / 'runs/L4/spin'
target = root / 'runs/L4/spin_rerun'
if target.exists():
    raise SystemExit(f'Refusing to overwrite {target}')
target.mkdir()
for name in ('gotm.yaml', 'fabm.yaml', 'runrun.sh', 'submit.slurm'):
    shutil.copy2(source/name, target/name)
for name in ('forcings', 'relax'):
    shutil.copytree(source/name, target/name, symlinks=False)
(target/'logs').mkdir()
print('Fresh deterministic spin directory:', target)
print('Run: cd', target)
print('Then: sbatch submit.slurm')
print('The 500-member assimilation continues to use the supplied spin/restart.nc.')

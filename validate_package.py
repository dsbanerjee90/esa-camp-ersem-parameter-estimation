# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).
from pathlib import Path
import ast
import subprocess
import yaml

root = Path(__file__).resolve().parent
experiment = root / 'runs/L4/assim/ens500_science_v3_T090_F010'
spin = root / 'runs/L4/spin'

def contained(path):
    try:
        path.resolve().relative_to(root)
    except ValueError:
        raise RuntimeError(f'External path dependency: {path}')

def require_file(base, value):
    path = Path(value)
    if path.is_absolute():
        raise RuntimeError(f'Absolute input path is not portable: {path}')
    path = base / path
    contained(path)
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f'Missing/empty input: {path}')

def walk_inputs(node, base):
    if isinstance(node, dict):
        # Reject absolute file references even when currently inactive.
        for key, value in node.items():
            if key in ('file', 'yaml_file') and isinstance(value, str) and value:
                if Path(value).is_absolute():
                    raise RuntimeError(f'Absolute file reference in {base}: {value}')
            if key == 'yaml_file' and value:
                require_file(base, value)
        if node.get('method') == 'file':
            require_file(base, node['file'])
        for value in node.values():
            walk_inputs(value, base)
    elif isinstance(node, list):
        for value in node:
            walk_inputs(value, base)

for path in root.rglob('*'):
    if path.is_symlink():
        raise RuntimeError(f'Symlink found: {path}')
    if not path.is_file():
        continue
    if path.suffix == '.py':
        ast.parse(path.read_text(), filename=str(path))
    if path.suffix in ('.sh', '.slurm'):
        subprocess.run(['bash', '-n', str(path)], check=True)
    if path.suffix in ('.sh', '.slurm', '.py') and 'runs' in path.relative_to(root).parts:
        text = path.read_text()
        if '/work/dba/' in text or '/scyllapfs/scratch/dba/' in text:
            raise RuntimeError(f'Original-user path in runtime code: {path}')

configs = []
for base in (spin, experiment):
    cfg = yaml.safe_load((base/'gotm.yaml').read_text())
    assert cfg['location']['name'] == 'L4'
    walk_inputs(cfg, base)
    walk_inputs(yaml.safe_load((base/'fabm.yaml').read_text()), base)
    configs.append(cfg)
assert configs[0]['time']['stop'] == configs[1]['time']['start'], 'Spin/assim times do not join'
assert configs[0]['restart']['load'] is False
assert configs[1]['restart']['load'] is True
for stream in yaml.safe_load((experiment/'config/obs.yaml').read_text())['streams']:
    require_file(experiment, stream['file'])
require_file(experiment, '../../spin/restart.nc')
settings = (experiment/'experiment.env').read_text()
assert 'export N_ENS=500\n' in settings and 'export MPI_GOTM=500\n' in settings
assert '#SBATCH --ntasks=501\n' in (experiment/'submit.slurm').read_text()
print('Portable package checks passed: no symlinks, internal inputs, matching dates, 500 members / 501 tasks.')

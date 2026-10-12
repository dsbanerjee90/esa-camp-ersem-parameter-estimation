# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA) under Climate and Marine Projection (CAMP) Project.
# If you use the code mention our Zenodo repo.

import ast
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import netCDF4
import eatpy
from parameter_prior import load_centres, prior_range
from parameter_evolution_plugin import _read_observation_schedule

base = Path(__file__).resolve().parent
import os
os.chdir(base)
for file in base.glob('*.py'):
    ast.parse(file.read_text(), filename=str(file))
gotm = yaml.safe_load(Path('gotm.yaml').read_text())
assert gotm['location']['name'] == 'L4', 'Not the L4 configuration'
assert gotm['restart']['load'] is True, 'Restart loading must be enabled'
start, stop = (pd.Timestamp(gotm['time'][k]) for k in ('start', 'stop'))
assert start < stop
centres = load_centres('fabm.yaml')
for (_, ptype), centre in centres.items():
    prior_range(centre, ptype)
print(f'L4: {start} to {stop}; depth={gotm["location"]["depth"]} m; levels={gotm["grid"]["nlev"]}')
print(f'Validated {len(centres)} parameter priors from L4 fabm.yaml')

cfg = yaml.safe_load(Path('config/obs.yaml').read_text())
allowed = {'O2_o', 'total_chlorophyll', 'N3_n', 'N1_p', 'N5_s', 'temp', 'salt'}
assert {s['var'] for s in cfg['streams']} == allowed
all_times = set()
for stream in cfg['streams']:
    path = Path(stream['file'])
    times = []
    rows = 0
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if len(fields) != 5:
            raise ValueError(f'{path}:{lineno}: expected date time depth value sd')
        t = pd.Timestamp(' '.join(fields[:2]))
        depth, value, sd = map(float, fields[2:])
        if not np.isfinite([depth, value, sd]).all() or sd <= 0:
            raise ValueError(f'{path}:{lineno}: invalid value or uncertainty')
        if stream['var'] == 'total_chlorophyll' and value <= 0:
            raise ValueError(f'{path}:{lineno}: chlorophyll must be positive for log10 DA')
        if not start <= t <= stop:
            raise ValueError(f'{path}:{lineno}: observation outside model period: {t}')
        if times and t < times[-1]:
            raise ValueError(f'{path}:{lineno}: timestamps are not sorted')
        times.append(t)
        all_times.add(t)
        rows += 1
    if not rows:
        raise ValueError(f'Empty observation stream: {path}')
    print(f'{stream["var"]}: {rows} rows; {len(set(times))} times; {times[0]} to {times[-1]}')
schedule = _read_observation_schedule('config/obs.yaml')
assert list(schedule) == sorted(all_times)
gaps = np.diff(schedule.asi8) / (86400 * 1e9)
print(f'Union schedule: {len(schedule)} times; gap range {gaps.min():.8g} to {gaps.max():.8g} days')
dt = float(gotm['time']['dt'])
off_grid = sum(not np.isclose(((t-start).total_seconds()/dt) % 1, 0, atol=1e-8)
               for t in schedule)
print(f'DEEP: Observation times off nominal {dt:g}-second grid: {off_grid}')
print('DEEP: Original timestamps are preserved. V3 requires callback times to match these timestamps exactly.')

def check_inputs(node):
    if isinstance(node, dict):
        if node.get('method') == 'file':
            p = Path(node['file'])
            if not p.is_file() or p.stat().st_size == 0:
                raise ValueError(f'Missing active input: {p}')
        for val in node.values():
            check_inputs(val)
    elif isinstance(node, list):
        for val in node:
            check_inputs(val)
check_inputs(gotm)
restart = Path('../../spin/restart.nc').resolve()
assert restart.is_file() and restart.stat().st_size > 0
with netCDF4.Dataset(restart) as nc:
    print('Spin restart:', restart)
    for name, var in nc.variables.items():
        data = np.ma.asarray(var[:])
        if np.ma.is_masked(data):
            raise ValueError(f'Masked restart data: {name}')
        values = np.asarray(data)
        if np.issubdtype(values.dtype, np.number):
            if not np.isfinite(values).all():
                raise ValueError(f'Non-finite restart data: {name}')
            if np.any(values < 0) and name not in {'h', 'z', 'zi', 'lon', 'lat'}:
                print(f'  Sign-preserving perturbation: {name}; min={values.min():.8g}')
        if name == 'time' and hasattr(var, 'units'):
            print('Restart time:', netCDF4.num2date(var[:], var.units,
                  calendar=getattr(var, 'calendar', 'standard')))
print('DEEP: Yay! Preflight passed. This is an input check though, not a completed MPI/model test.')

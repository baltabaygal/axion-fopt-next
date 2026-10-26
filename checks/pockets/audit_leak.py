"""Read-only snapshot audit of existing leak runs, including unfinished runs.

Writes only to a fresh output directory. No absent sample is treated as zero,
and far-field storage is never interpreted as cumulative escaped energy.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from model import read_slice
from run_prototype import wall_cells


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(source):
    contract = json.loads((source/'contract.json').read_text())
    reports = {}
    for name, spec in contract['runs'].items():
        directory = source/name
        path = directory/'native/history.csv'
        if not path.exists():
            reports[name] = {'status': 'missing_history'}
            continue
        raw = path.read_bytes()
        lines = raw.splitlines(keepends=True)
        columns = len(lines[0].strip().split(b','))
        # A running process may be in the middle of writing its final row.
        complete = [line for line in lines[1:] if line.endswith(b'\n')
                    and len(line.strip().split(b',')) == columns]
        if not complete:
            reports[name] = {'status': 'no_complete_samples'}
            continue
        h = np.atleast_1d(np.genfromtxt(io.BytesIO(lines[0]+b''.join(complete)), delimiter=',', names=True))
        if not all(np.isfinite(h[k]).all() for k in h.dtype.names):
            raise ValueError(f'Nonfinite history: {path}')
        t, r = h['u'], h['fv_radius']
        dead = np.flatnonzero(r < contract['gates']['survival_radius'])
        life = float(t[dead[0]]) if len(dead) else None
        alive = t < life if life is not None else np.ones(len(t), bool)
        ax0 = sum(float(h[k][0]) for k in ('axion_kinetic', 'axion_gradient', 'coupling'))
        r0 = contract['r0']; box = contract['box']
        report = dict(status='completed' if (directory/'native/final.bin').exists() else 'partial',
                      last_time=float(t[-1]), final_radius=float(r[-1]),
                      lifetime=life, censored=life is None,
                      energy_drift_alive=float(np.max(np.abs(h['total'][alive]/h['total'][0]-1))),
                      history_sha256=hashlib.sha256(raw).hexdigest(),
                      requested_duration=spec['duration'],
                      earliest_causal_image_return=box-2*r0,
                      note='Image-return bound uses initial radius and speed <=1; it is not a measured contamination time.')
        widths = []
        unmeasurable_widths = 0
        for p in sorted((directory/'native').glob('slice_*.bin')):
            meta, s = read_slice(p)
            if meta['u'] <= t[-1] and (life is None or meta['u'] < life):
                w = wall_cells(s, box/spec['n'])
                if np.isfinite(w):
                    widths.append(float(w))
                else:
                    unmeasurable_widths += 1
        report['minimum_wall_cells_alive'] = min(widths) if widths else None
        report['unmeasurable_wall_slices_alive'] = unmeasurable_widths
        report['energy_gate'] = report['energy_drift_alive'] < contract['gates']['energy_drift']
        if ax0 > 0:
            valid = h['fv_volume'] > 0
            far = h['axion_far']/ax0
            report['far_fraction_at_100'] = (float(far[np.flatnonzero(np.isclose(t, 100., atol=1e-8, rtol=0))[0]])
                                            if np.any(np.isclose(t, 100., atol=1e-8, rtol=0) & valid) else None)
            lo, hi = contract['classes']['leak_window']
            hi = min(hi, float(t[-1]), life if life is not None else np.inf)
            mask = (t >= lo) & (t <= hi) & valid
            sufficient = mask.sum() >= 2 and np.ptp(t[mask]) >= contract['classes']['minimum_window']
            report['far_storage_slope'] = float(np.polyfit(t[mask], far[mask], 1)[0]) if sufficient else None
            report['slope_window'] = [lo, hi]
            report['initial_far_fraction'] = float(far[0])
            report['slope_causally_isolated'] = hi < box-2*r0
        reports[name] = report
    files = contract.get('sources', {})
    matches = {name: digest(Path(__file__).parent/name) == sha for name, sha in files.items()}
    return dict(source=str(source.resolve()), contract_sha256=digest(source/'contract.json'),
                current_sources_match=matches, runs=reports,
                interpretation='Snapshot audit only. Partial runs retain observed threshold crossings; '
                'survival beyond the last saved sample is unknown. '
                'Far storage slopes are not irreversible loss rates. No continuum or isolated-pocket claim.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be fresh')
    result = audit(args.source)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()

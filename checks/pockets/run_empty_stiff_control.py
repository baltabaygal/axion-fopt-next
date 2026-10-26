"""Matched compiled empty control; EMPTY_STIFF_CONTROL.md."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import config,eigen_state,history,radial_collapse,read_slice,read_state,run
from wall_diagnostics import local_widths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    exe=Path('build/pockets-omp/pocket_evolve').resolve()
    validation=json.loads(Path('outputs/pockets/validation-omp-02/validation.json').read_text())
    if not validation['passed'] or hashlib.sha256(exe.read_bytes()).hexdigest()!=validation['executable_sha256']:
        raise ValueError('Validated accurate OpenMP executable required')
    s=eigen_state(256,64.,20.,1.)
    _,gas_initial=read_state(Path('outputs/pockets/stiff-wall-01/n256/native/initial.bin'))
    matched=bool(np.array_equal(s[...,0],gas_initial[...,0]) and np.array_equal(s[...,2],gas_initial[...,2]))
    del gas_initial
    cfg=config(256,64.,.0125,20.,lam=8.,dv=.01,mass_tv=2.,mass_power=4,every=20,slices=80,frames=800)
    parent=Path(__file__).parent
    (args.output/'contract.json').write_text(json.dumps(dict(config=cfg,exact_matched_dilaton=matched,
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        sources={f:hashlib.sha256((parent/f).read_bytes()).hexdigest() for f in (
            'EMPTY_STIFF_CONTROL.md','run_empty_stiff_control.py','model.py','wall_diagnostics.py')}),indent=2)+'\n')
    error=None
    try:
        run(s,cfg,args.output/'empty',exe)
    except RuntimeError as exc:
        error=str(exc)
    h=history(args.output/'empty');_,initial=read_state(args.output/'empty/native/initial.bin')
    transfer=bool(np.array_equal(s,initial));del s,initial
    t,r=radial_collapse(20.,1.,8.,.01,32.,.025,20.,.002,125)
    np.savetxt(args.output/'radial.csv',np.c_[t,r],delimiter=',',header='u,radius',comments='')
    first=np.flatnonzero(r<15);trusted_end=t[first[0]-1] if len(first) else t[-1]
    mask=h['u']<=trusted_end
    difference=float(np.max(np.abs(h['fv_radius'][mask]-np.interp(h['u'][mask],t,r))))
    drift=float(np.max(np.abs(h['total'][mask]/h['total'][0]-1)))
    widths=[]
    for path in sorted((args.output/'empty/native').glob('slice_*.bin')):
        meta,plane=read_slice(path)
        if meta['u']<=trusted_end:
            widths.append(dict(time=meta['u'],**local_widths(plane[...,0],.25)))
    gas=history(Path('outputs/pockets/stiff-wall-01/n256'))
    crossing=float(t[first[0]]) if len(first) else None
    gates=dict(matched=matched,transfer=transfer,reference=difference<.25,energy=drift<5e-4,
               trusted_completed=bool(h['u'][-1]>=trusted_end),wall=bool(widths) and all(w['diagnostic_pass'] for w in widths))
    report=dict(gates=gates,solver_error=error,trusted_end=float(trusted_end),
        maximum_radius_error=difference,energy_drift_trusted=drift,widths_trusted=widths,
        reference_R15_crossing=crossing,gas_radius_at_crossing=float(np.interp(crossing,gas['u'],gas['fv_radius'])) if crossing else None,
        final_recorded_radius=float(h['fv_radius'][-1]),last_time=float(h['u'][-1]),
        interpretation='Early controlled contraction only; unresolved later implosion is not an accurate collapse-time result.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='widths_trusted'},indent=2))


if __name__=='__main__':
    main()

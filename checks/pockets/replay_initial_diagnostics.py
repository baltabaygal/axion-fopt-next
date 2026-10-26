"""Fresh one-step replay to verify compensated initial-energy diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from model import energies,history,read_state


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='Completed original case directory')
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((args.source/'input.json').read_text())
    original_input=Path(cfg['initial_state'])
    cfg.update(tMax=cfg['dt'],output_every=1,snapshot_every=0,slice_every=0,
               result_dir=str((args.output/'native').resolve()))
    (args.output/'input.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (args.output/'input.in').write_text('\n'.join(
        f'{k} = {str(v).lower() if isinstance(v,bool) else v}' for k,v in cfg.items())+'\n')
    parent=Path(__file__).parent
    provenance=dict(source=str(args.source.resolve()),source_input_sha256=sha(original_input),
                     executable=str(args.executable.resolve()),executable_sha256=sha(args.executable),
                     sources={f:sha(parent/f) for f in ('DIAGNOSTIC_PRECISION.md','accurate_sum.h',
                               'evolve.cpp','pocket_monitor.h','model.py','replay_initial_diagnostics.py')})
    (args.output/'contract.json').write_text(json.dumps(provenance,indent=2)+'\n')
    result=subprocess.run([str(args.executable.resolve()),'input='+str((args.output/'input.in').resolve())],
                          capture_output=True,text=True)
    (args.output/'run.log').write_text(result.stdout+result.stderr)
    result.check_returncode()
    meta,initial=read_state(args.output/'native/initial.bin')
    expected=np.memmap(original_input,dtype='<f8',mode='r',offset=20,shape=initial.shape)
    transfer=bool(np.array_equal(initial,expected));del expected
    independent=energies(initial,meta['box'],cfg['lambda'],cfg['delta_v'],cfg['mass_tv']**2,
                         far_margin=cfg.get('far_margin',6.),mass_power=cfg.get('mass_power',2))
    updated=history(args.output);old=history(args.source);scale=abs(independent['total'])
    energy_keys=('total','dilaton_kinetic','dilaton_gradient','dilaton_potential',
                 'axion_kinetic','axion_gradient','coupling','axion_fv','axion_wall','axion_tv','axion_far')
    errors={k:float(abs(updated[k][0]-independent[k])/scale) for k in energy_keys}
    old_errors={k:float(abs(old[k][0]-independent[k])/scale) for k in energy_keys}
    axion=sum(independent[k] for k in ('axion_kinetic','axion_gradient','coupling'))
    new_axion=sum(float(updated[k][0]) for k in ('axion_kinetic','axion_gradient','coupling'))
    old_axion=sum(float(old[k][0]) for k in ('axion_kinetic','axion_gradient','coupling'))
    report=dict(exact_transfer=transfer,max_initial_energy_error=max(errors.values()),
                old_max_initial_energy_error=max(old_errors.values()),component_errors=errors,
                original_component_errors=old_errors,
                axion_relative_error=float(abs(new_axion-axion)/axion),
                old_axion_relative_error=float(abs(old_axion-axion)/axion),
                passed=bool(transfer and max(errors.values())<1e-12 and abs(new_axion-axion)/axion<1e-12),
                interpretation='Initial-state diagnostic verification only; original failed report retained.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()

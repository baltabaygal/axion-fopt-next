"""Static checkpoint and cross-backend verification; RESTART_CONTRACT.md."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model import config,read_state,run
from validate import small_state


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--serial',type=Path,default=Path('build/pockets-serial-check/pocket_evolve'))
    p.add_argument('--threaded',type=Path,default=Path('build/pockets-omp/pocket_evolve'))
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    parent=Path(__file__).parent
    (args.output/'contract.json').write_text(json.dumps(dict(
        serial_sha256=sha(args.serial),threaded_sha256=sha(args.threaded),
        sources={f:sha(parent/f) for f in ('RESTART_CONTRACT.md','verify_restart.py','model.py','validate.py')},
        relative_state_tolerance=1e-12),indent=2)+'\n')
    s=small_state(6,6.)
    cfg=lambda duration:config(6,6.,.005,duration,lam=8.,dv=.01,mass_tv=2.,mass_power=4,every=20)
    run(s,cfg(1.),args.output/'full',args.serial.resolve())
    run(s,cfg(.4),args.output/'prefix',args.serial.resolve())
    _,checkpoint=read_state(args.output/'prefix/native/final.bin')
    run(checkpoint,cfg(.6),args.output/'continuation',args.threaded.resolve())
    _,transferred=read_state(args.output/'continuation/native/initial.bin')
    _,full=read_state(args.output/'full/native/final.bin')
    _,continued=read_state(args.output/'continuation/native/final.bin')
    error=float(np.max(np.abs(full-continued))/max(1.,np.abs(full).max()))
    transfer=bool(np.array_equal(checkpoint,transferred))
    report=dict(exact_checkpoint_transfer=transfer,relative_final_state_error=error,
                bitwise_final_agreement=bool(np.array_equal(full,continued)),
                passed=bool(transfer and error<1e-12),
                interpretation='Static canonical state restart only; no physical survival or expanding restart claim.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()

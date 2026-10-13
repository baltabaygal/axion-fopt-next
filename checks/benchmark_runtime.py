"""Alternating before/after 128^3 runs, with unchanged physics and outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"checks/single_bubble"))
from reference import read_state, observables


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline",type=Path,default=ROOT/"outputs/runtime-128-baseline/profiled_evolve")
    p.add_argument("--optimized",type=Path,default=ROOT/"build/cosmolattice/axion_fopt_evolve")
    p.add_argument("--output",type=Path,default=ROOT/"outputs/runtime-128-optimized")
    p.add_argument("--repeats",type=int,default=3)
    args=p.parse_args()
    if args.repeats<1: p.error("repeats must be positive")
    executables={"baseline":args.baseline.resolve(),"optimized":args.optimized.resolve()}
    for path in executables.values():
        if not path.is_file(): p.error(f"Missing executable {path}")
    root=args.output.resolve(); root.mkdir(parents=True,exist_ok=False)
    original=ROOT/"outputs/lee-single-bubble/N128_v0.5_ref1"
    _,golden=read_state(original/"final.bin")
    golden_history=np.genfromtxt(original/"history.csv",delimiter=",",names=True)
    golden_observables,_=observables(golden[...,0],golden[...,1],2*np.pi)
    params=json.loads((ROOT/"outputs/lee-single-bubble/N128_v0.5_ref1.input.json").read_text())
    params.pop("snapshot_at",None)
    params["profile"]="true"
    records=[]
    for repeat in range(args.repeats):
        for label in (("baseline","optimized") if repeat%2==0 else ("optimized","baseline")):
            out=root/f"{label}_{repeat+1}"
            run_params=params | {"result_dir":str(out)}
            (root/f"{out.name}.input.json").write_text(json.dumps(run_params,indent=2)+"\n")
            print(f"START {out.name}",flush=True)
            begin=time.monotonic()
            command=[str(executables[label]),*[f"{k}={v}" for k,v in run_params.items()]]
            # wait4 captures each individual child's RSS without sysctl-based
            # platform timing tools. FFT/comparison below is outside this timer.
            with (root/f"{out.name}.log").open("w") as log:
                proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
                _,exit_status,usage=os.wait4(proc.pid,0)
                proc.returncode=os.waitstatus_to_exitcode(exit_status)
            elapsed=time.monotonic()-begin
            if proc.returncode: raise RuntimeError(f"Failed {out.name}; see log")
            _,state=read_state(out/"final.bin")
            history=np.genfromtxt(out/"history.csv",delimiter=",",names=True)
            np.testing.assert_allclose(state,golden,atol=1e-12,rtol=0)
            history_error=0.
            for column in history.dtype.names:
                delta=float(np.max(np.abs(history[column]-golden_history[column])))
                history_error=max(history_error,delta)
                np.testing.assert_allclose(history[column],golden_history[column],atol=1e-11,rtol=0)
            stats,_=observables(state[...,0],state[...,1],2*np.pi)
            number_error=abs(stats["number_continuum"]/golden_observables["number_continuum"]-1)
            if number_error>=1e-11: raise RuntimeError("Particle number changed")
            record=dict(label=label,repeat=repeat+1,wall_seconds=elapsed,
                        final_snapshot_sha256=hashlib.sha256((out/"final.bin").read_bytes()).hexdigest(),
                        peak_rss_bytes=usage.ru_maxrss if sys.platform=="darwin" else usage.ru_maxrss*1024,
                        timings=json.loads((out/"timings.json").read_text()) if (out/"timings.json").exists() else None,
                        final_state_max_error=float(np.max(np.abs(state-golden))),
                        history_max_error=history_error,number_relative_error=float(number_error),
                        number_continuum=stats["number_continuum"])
            records.append(record)
            (root/"runs.json").write_text(json.dumps(records,indent=2)+"\n")
            print(f"DONE {out.name}: {elapsed:.3f}s; final state error {record['final_state_max_error']:.3e}",flush=True)
    medians={label:float(np.median([r["wall_seconds"] for r in records if r["label"]==label])) for label in executables}
    source_files=[Path(__file__),*[ROOT/"adapters/cosmolattice"/name for name in
                  ("axion_fopt.h","potential_functions.h","evolve.cpp","prepare_overlay.py","CMakeLists.txt")]]
    result=dict(status="passed",configuration=params,repeats=args.repeats,runs=records,
                median_wall_seconds=medians,speedup=medians["baseline"]/medians["optimized"],
                executable_sha256={label:hashlib.sha256(path.read_bytes()).hexdigest() for label,path in executables.items()},
                source_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
                golden_snapshot_sha256=hashlib.sha256((original/"final.bin").read_bytes()).hexdigest(),
                methodology="3 alternating pairs by default; full executable wall time including diagnostics and two snapshots; FFT and comparisons excluded; double precision; identical 1280 steps and output_every=2",
                timing_note="mass_update_seconds overlaps initialization/evolution/diagnostics; do not add it to those phase totals")
    (root/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"medians":medians,"speedup":result["speedup"]},indent=2))


if __name__=="__main__": main()

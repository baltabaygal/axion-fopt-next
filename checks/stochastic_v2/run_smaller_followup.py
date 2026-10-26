"""One prospectively specified L20 follow-up; no parameter scan or damping."""
import argparse
import json
from pathlib import Path
import shutil
from model import seed,run,PARAMS
from run_small_box import analyze
from validate_backend import verify_build
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();build=verify_build(args.executable)
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(Path(__file__).with_name('SMALL_BOX_FOLLOWUP.md'),out/'contract.md')
    (out/'parameters.json').write_text(json.dumps(dict(box=20.,n=80,dx=.25,dt=.01,end=32.,physical=PARAMS,R0=8.,w=2.,theta0=1.,initial_velocities=[0,0],friction=0,expansion=False),indent=2)+'\n')
    seconds=run(args.executable,out/'main',seed(80,20.,1.),20.,.01,32.,cadence=.04,frames=.2,timeout=180)
    r=analyze(out/'main',80,.01,box=20.,end=32.)
    report=dict(results=r,passed=bool(all(r['checks'].values())),goal_reached=r['first_sustained_full_interval'] is not None,
        seconds=seconds,build=build,source_sha256={str(f.resolve()):sha(f) for f in [Path(__file__),Path(__file__).with_name('run_small_box.py'),Path(__file__).with_name('model.py')]})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='build'},indent=2))

if __name__=='__main__':main()

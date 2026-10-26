"""Build only our owned driver; reuse existing serial CosmoLattice libraries."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--base-manifest',type=Path,default=ROOT/'build/pockets-formation-light-01/build.json')
    args=p.parse_args();dest=args.output.resolve();dest.mkdir(parents=True,exist_ok=False)
    base=json.loads(args.base_manifest.read_text());compile=list(base['compile']);link=list(base['link'])
    old_obj=compile[compile.index('-o')+1]
    compile[compile.index('-c')+1]=str(HERE/'stochastic_evolve.cpp');compile[compile.index('-o')+1]=str(dest/'evolve.o')
    compile=['-O2' if x.startswith('-O') else x for x in compile]
    link=[str(dest/'evolve.o') if x==old_obj else ('-O2' if x.startswith('-O') else x) for x in link]
    link[link.index('-o')+1]=str(dest/'bubble_evolve')
    library_build=Path(base['existing_library_build']);started=time.monotonic()
    for name,cmd,cwd in [('compile',compile,HERE),('link',link,library_build)]:
        result=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,timeout=45)
        (dest/(name+'.log')).write_text(result.stdout+result.stderr)
        result.check_returncode()
    deps=[HERE/n for n in ['stochastic_evolve.cpp','bubble_model.h','state_io.h','potential.h','bubble_insertion.h','build_stochastic_backend.py']]
    deps += [HERE/'../pockets/initial_state.h',HERE/'../pockets/accurate_sum.h']
    libraries=[Path(x) if Path(x).is_absolute() else library_build/x for x in link if x.endswith(('.a','.dylib'))]
    metadata=dict(compile=compile,link=link,library_build=str(library_build),
        source_sha256={str(x.resolve()):sha(x) for x in deps},library_sha256={str(x):sha(x) for x in libraries},
        base_manifest_sha256=sha(args.base_manifest),executable_sha256=sha(dest/'bubble_evolve'),
        seconds=time.monotonic()-started,scope='3D coupled fields; prescribed RD expansion; actual CosmoLattice PV2, one CPU thread')
    (dest/'build.json').write_text(json.dumps(metadata,indent=2)+'\n');print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()

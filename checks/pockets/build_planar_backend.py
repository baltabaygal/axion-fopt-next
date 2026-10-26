"""Compile only the planar driver against existing CPU libraries, serially."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-manifest',type=Path,required=True)
    parser.add_argument('--library-build',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--optimization',choices=('O0','O2'),default='O2')
    args=parser.parse_args();base=json.loads(args.base_manifest.read_text());dest=args.output.resolve();dest.mkdir(parents=True,exist_ok=False)
    for name,digest in base['source_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('Base source changed: '+name)
    compile=list(base['compile']);old_object=compile[compile.index('-o')+1]
    compile[compile.index('-o')+1]=str(dest/'line.o')
    compile=[('-'+args.optimization) if arg in ('-O0','-O1','-O2','-O3') else arg for arg in compile]
    link=[str(dest/'line.o') if arg==old_object else arg for arg in base['link']]
    link[link.index('-o')+1]=str(dest/'pocket_line_evolve')
    link=[('-'+args.optimization) if arg in ('-O0','-O1','-O2','-O3') else arg for arg in link]
    started=time.monotonic()
    for label,command,cwd in (('compile',compile,HERE),('link',link,args.library_build.resolve())):
        result=subprocess.run(command,cwd=cwd,capture_output=True,text=True,timeout=45)
        (dest/(label+'.log')).write_text(result.stdout+result.stderr);result.check_returncode()
    libraries={}
    for arg in link:
        if arg.endswith(('.a','.dylib')):
            path=Path(arg);path=path if path.is_absolute() else args.library_build.resolve()/path
            libraries[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    report=dict(compile=compile,link=link,library_build=str(args.library_build.resolve()),
        library_sha256=libraries,optimization=args.optimization,fast_math_added=False,
        elapsed_seconds=time.monotonic()-started,source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in base['source_sha256']},
        executable_sha256=hashlib.sha256((dest/'pocket_line_evolve').read_bytes()).hexdigest(),
        base_manifest_sha256=hashlib.sha256(args.base_manifest.read_bytes()).hexdigest(),
        build_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Actual NDim=1 coupled fields; optional diagnostics, no field-model change')
    (dest/'build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(optimization=report['optimization'],elapsed_seconds=report['elapsed_seconds'],executable=str(dest/'pocket_line_evolve')),indent=2))

if __name__=='__main__':main()

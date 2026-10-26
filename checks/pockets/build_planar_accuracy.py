"""Serial build of the optional planar stencil hook, against existing libraries."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from prepare_line_accuracy_overlay import prepare_overlay

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-manifest',type=Path,required=True)
    p.add_argument('--library-build',type=Path,required=True)
    p.add_argument('--vendor',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();base=json.loads(a.base_manifest.read_text());dest=a.output.resolve();dest.mkdir(parents=True,exist_ok=False)
    commit=subprocess.run(['git','-C',str(a.vendor),'rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()
    if commit!='0dad16ec2c894a7c6abbea0e98734126b3470b36':raise RuntimeError('Unsupported CosmoLattice commit')
    overlay=prepare_overlay(a.vendor,dest/'overlay')
    compile=list(base['compile']);old_object=compile[compile.index('-o')+1];compile[compile.index('-o')+1]=str(dest/'line.o')
    compile=[arg if arg not in ('-O0','-O1','-O3') else '-O2' for arg in compile]
    compile[1:1]=['-DPOCKETS_LINE_CUSTOM_LAPLACIAN=1','-I'+str(dest/'overlay'),'-I'+str(HERE)]
    link=[str(dest/'line.o') if arg==old_object else arg for arg in base['link']]
    link[link.index('-o')+1]=str(dest/'pocket_line_evolve')
    started=time.monotonic()
    for label,command,cwd in (('compile',compile,HERE),('link',link,a.library_build.resolve())):
        result=subprocess.run(command,cwd=cwd,capture_output=True,text=True,timeout=45)
        (dest/(label+'.log')).write_text(result.stdout+result.stderr)
        if result.returncode:print(result.stderr[-5000:]);result.check_returncode()
    libraries={}
    for arg in link:
        if arg.endswith(('.a','.dylib')):
            path=Path(arg);path=path if path.is_absolute() else a.library_build.resolve()/path
            libraries[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    names=set(base['source_sha256'])|{'line_laplacian.h','line_spatial_energy.h','line_wall_geometry.h','prepare_line_accuracy_overlay.py'}
    sources={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sorted(names)}
    sources[overlay['overlay']]=overlay['overlay_sha256']
    report=dict(compile=compile,link=link,library_build=str(a.library_build.resolve()),library_sha256=libraries,
        optimization='O2',fast_math_added=False,elapsed_seconds=time.monotonic()-started,
        source_sha256=sources,overlay=overlay,cosmolattice_commit=commit,
        supports_line_spatial_orders=[2,4],supports_line_wall_geometry=True,ghost_cells={'2':1,'4':2},
        executable_sha256=hashlib.sha256((dest/'pocket_line_evolve').read_bytes()).hexdigest(),
        base_manifest_sha256=hashlib.sha256(a.base_manifest.read_bytes()).hexdigest(),
        build_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Same canonical planar fields and PV2; optional consistent order-4 spatial Hamiltonian')
    (dest/'build.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(orders=report['supports_line_spatial_orders'],elapsed_seconds=report['elapsed_seconds']),indent=2))

if __name__=='__main__':main()

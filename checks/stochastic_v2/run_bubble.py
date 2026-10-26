"""Bounded 3D single-bubble pilot and matched numerical/empty controls."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from scipy.optimize import brentq
from model import PARAMS,seed,read_state,run,terms
from build_backend import sha
from validate_backend import verify_build

CASES=[('N48-dt004',48,48.,.04,1.),('N48-dt002',48,48.,.02,1.),
       ('N48-dt001',48,48.,.01,1.),('N64',64,48.,.02,1.),
       ('N96',96,48.,.02,1.),('empty-N64',64,48.,.02,0.),('boundary-N80',80,60.,.02,1.)]

def contour(line,box,transverse2):
    n=len(line);x=(np.arange(n)+.5)*box/n-box/2
    r=np.sqrt(x[n//2:]**2+transverse2);s=line[n//2:]
    result={};counts={}
    for value,name in ((.1,'r10'),(.5,'radius'),(.9,'r90')):
        indices=np.nonzero((s[:-1]>=value)&(s[1:]<value))[0];counts[name]=len(indices)
        if len(indices)!=1:result[name]=None
        else:
            i=int(indices[0]);result[name]=float(r[i]+(r[i+1]-r[i])*(value-s[i])/(s[i+1]-s[i]))
    result['crossing_counts']=counts
    result['width']=None if result['r10'] is None or result['r90'] is None else result['r10']-result['r90']
    return result

def axis_line(field,axis):
    # Average the four nearest lines to a central coordinate axis. For a radial
    # field all have r^2=x^2+dx^2/2; averaging also diagnoses axis anisotropy.
    n=len(field);m=np.moveaxis(field,axis,0)
    return m[:,n//2-1:n//2+1,n//2-1:n//2+1].mean(axis=(1,2))

def analyze(directory,n,box,dt,theta):
    native=directory/'native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    _,initial=read_state(native/'initial.bin');_,final=read_state(native/'final.bin');dx=box/n
    energy={k:float(v.sum()*dx**3) for k,v in terms(final,dx).items()}
    snap_error=max(abs(energy[k]-float(h[k][-1])) for k in energy)/abs(float(h['total'][0]))
    background=PARAMS['dv']*box**3;excess=float(h['total'][0]-background)
    delta=np.max(np.abs(h['total']-h['total'][0]))
    q=np.r_[0,np.cumsum((h['Q_A'][1:]+h['Q_A'][:-1])/2*np.diff(h['t']))]
    exchange=float(np.max(np.abs(h['E_A']-h['E_A'][0]-q))/max(float(h['E_A'][0]),1e-30))
    trajectory=[]
    for path in sorted(native.glob('frame_*.bin')):
        meta,plane=read_state(path)
        line=plane[:,n//2-1:n//2+1,0].mean(axis=1)
        c=contour(line,box,dx*dx/2);c['t']=meta['t'];trajectory.append(c)
    endpoint=[contour(axis_line(final[...,0],i),box,dx*dx/2) for i in range(3)]
    radii=[c['radius'] for c in endpoint]
    return dict(n=n,box=box,dt=dt,theta_initial=theta,initial_transfer_exact=bool(np.array_equal(initial,seed(n,box,theta))),
        initial_energy=float(h['total'][0]),initial_excess_above_fv=excess,
        initial_axion_energy=float(h['E_A'][0]),final_axion_energy=float(h['E_A'][-1]),
        final_radius=radii[0],radius_anisotropy=float(max(radii)-min(radii)) if None not in radii else None,
        final_axion_rms=float(np.sqrt(np.mean(final[...,1]**2))/PARAMS['fa']),
        energy_drift_total=float(delta/abs(h['total'][0])),energy_drift_excess=float(delta/max(abs(excess),1e-30)),
        snapshot_energy_relative=float(snap_error),exchange_residual=exchange,
        empty_exact=bool(np.all(final[...,1::2]==0)) if theta==0 else None,
        wall_trajectory=trajectory,final_axis_contours=endpoint)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--executable',type=Path,required=True)
    parser.add_argument('--validation',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--analyze-only',action='store_true');args=parser.parse_args()
    out=args.output.resolve();build=verify_build(args.executable)
    validation=json.loads(args.validation.read_text())
    if not validation['passed'] or validation['build']['executable_sha256']!=build['executable_sha256']:raise RuntimeError('Validated matching backend required')
    if not args.analyze_only:
        out.mkdir(parents=True,exist_ok=False)
        (out/'contract.md').write_text(Path(__file__).with_name('EXPERIMENT.md').read_text())
        (out/'build.json').write_text(json.dumps(build,indent=2)+'\n')
        sources=out/'sources';sources.mkdir()
        for source in Path(__file__).parent.iterdir():
            if source.suffix in ('.py','.h','.cpp','.md'):shutil.copy2(source,sources/source.name)
        completed=[]
        for name,n,box,dt,theta in CASES:
            seconds=run(args.executable,out/name,seed(n,box,theta),box,dt,6.)
            completed.append(dict(name=name,seconds=seconds));(out/'timing.json').write_text(json.dumps(completed,indent=2)+'\n')
            print(json.dumps(completed[-1]),flush=True)
    results={name:analyze(out/name,n,box,dt,theta) for name,n,box,dt,theta in CASES}
    finals=[read_state(out/name/'native/final.bin')[1] for name in ('N48-dt004','N48-dt002','N48-dt001')]
    scale=np.array([1.,PARAMS['fa'],1.,PARAMS['fa']]);diff=[float(np.max(abs(finals[i]-finals[i+1])/scale)) for i in (0,1)]
    a,b=results['N64'],results['N96']
    spatial_radius=abs(a['final_radius']-b['final_radius'])/b['final_radius']
    spatial_energy=abs(a['final_axion_energy']-b['final_axion_energy'])/b['final_axion_energy']
    _,small=read_state(out/'N64/native/final.bin');_,large=read_state(out/'boundary-N80/native/final.bin')
    # Same dx, grid aligned after removing 8 cells on each side of the larger box.
    aligned=large[8:-8,8:-8,8:-8];x=(np.arange(64)+.5)*.75-24
    mask=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2<15**2
    boundary_rms=float(np.sqrt(np.mean(((small[mask]-aligned[mask])/scale)**2)))
    boundary_radius=abs(a['final_radius']-results['boundary-N80']['final_radius'])/a['final_radius']
    widths=[v['width']/(48/96) for v in b['wall_trajectory'] if v['width'] is not None]
    checks=dict(backend_reference=validation['passed'],initial_transfer=all(v['initial_transfer_exact'] for v in results.values()),
        snapshot_energies=all(v['snapshot_energy_relative']<1e-11 for v in results.values()),
        total_energy=all(v['energy_drift_excess']<.005 for v in results.values()),
        axion_exchange=all(v['exchange_residual']<.01 for v in results.values() if v['theta_initial']!=0),
        empty_exact=results['empty-N64']['empty_exact'],time_refinement=diff[0]/diff[1]>3.5,
        spatial_radius=spatial_radius<.05,spatial_axion_energy=spatial_energy<.1,
        resolved_wall=len(widths)==len(b['wall_trajectory']) and min(widths)>=5,
        boundary_fields=boundary_rms<.001,boundary_radius=boundary_radius<.005)
    # Full-potential true minimum at the imposed initial angle, not U alone.
    true_min=brentq(lambda s:PARAMS['lam']/2*s*(1-s)*(1-2*s)-12*PARAMS['dv']*s*s*(1-s)+2*PARAMS['fa']**2*PARAMS['mass']**2*s*(1-np.cos(1)),.8,1)
    report=dict(checks=checks,passed=bool(all(checks.values())),results=results,
        growth_demonstrated=b['final_radius']-b['wall_trajectory'][0]['radius']>.1,
        time_max_scaled_differences=diff,time_difference_ratio=diff[0]/diff[1],
        spatial_radius_relative=spatial_radius,spatial_axion_energy_relative=spatial_energy,
        boundary_rms=boundary_rms,boundary_radius_relative=boundary_radius,min_finest_wall_cells=min(widths) if widths else None,
        full_true_minimum_at_theta1=true_min,validation_report_sha256=sha(args.validation),build=build,
        source_sha256={str(p.resolve()):sha(p) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.h','.cpp','.md')},
        scope='Seeded finite-time scalar-only 3D pilot; no terminal speed, tunnelling, cosmology, relic density or continuum limit established')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('results','build','source_sha256')},indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()

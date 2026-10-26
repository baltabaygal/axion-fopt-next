"""Tiny C++ monitor/reference comparison; no framework evolution."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
from model import radius,laplacian,force as polynomial_force
from paper_reference import Parameters,force as formation_force
from monitor_reference import monitor
from tracking_monitor_reference import tracking_monitor

HERE=Path(__file__).resolve().parent
NAMES=('energy','phi_energy','axion_energy','phi_flux','axion_flux','exchange')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    binary=(args.output/'monitor_probe').resolve()
    compiler=shutil.which('c++')
    if not compiler:raise RuntimeError('C++ compiler unavailable')
    command=[compiler,'-std=c++17','-O0',str(HERE/'formation_monitor_probe.cpp'),'-o',str(binary)]
    build=subprocess.run(command,capture_output=True,text=True,timeout=30)
    (args.output/'build.log').write_text(build.stdout+build.stderr);build.check_returncode()
    n=12;box=8.;r=radius(n,box);state=np.zeros((n,n,n,4))
    state[...,0]=.5*(1+np.tanh((r-1.7)/.6))
    state[...,1]=2.*np.exp(-r*r/3.)
    state[...,2]=-.1*(np.roll(state[...,0],-1,0)-np.roll(state[...,0],1,0))+.02*np.exp(-r*r/3.)
    state[...,3]=.15*np.exp(-r*r/3.)
    np.save(args.output/'initial.npy',state)
    gates=dict(reference_scaled_error=1e-12,continuity_scaled_error=1e-7)
    def compiled(s,mode,label):
        text=' '.join(format(v,'.17g') for v in s.ravel())
        out=subprocess.run([str(binary),str(n),str(mode)],input=text,text=True,capture_output=True,check=True,timeout=10)
        (args.output/(label+'.txt')).write_text(out.stdout)
        lines=out.stdout.splitlines()
        values=[np.fromstring(line,sep=' ') for line in lines]
        if len(values)!=2 or values[0].size!=6 or values[1].size!=21:raise RuntimeError('Bad monitor output')
        if not all(np.isfinite(v).all() for v in values):raise RuntimeError('Nonfinite monitor output')
        return values
    results=[]
    for mode,label in ((0,'formation'),(1,'polynomial2'),(2,'polynomial4-stiff')):
        p=Parameters();power=4 if mode==2 else 2;b=12. if mode==2 else 0.
        formation=p if mode==0 else None
        references=[monitor(state,box,2.7,8.,.01,1.,power,b,formation),
                    tracking_monitor(state,box,2.4,3.5,8.,.01,1.,power,b,formation)]
        reference_arrays=[np.array([references[0]['monitor_'+name] for name in NAMES]),
            np.array([references[1]['tracking_'+name] for name in NAMES]+
            [references[1][key] for key in ('tracking_lab_phi_flux','tracking_lab_axion_flux','tracking_phi_advection','tracking_axion_advection')]+
            references[1]['centre']+references[1]['velocity']+[references[1]['phase_volume']]+
            references[1]['centroid_coherence'])]
        values=compiled(state,mode,label)
        fp,fa=formation_force(state[...,0],state[...,1],p) if mode==0 else polynomial_force(state[...,0],state[...,1],8.,.01,1.,power,b)
        rhs=np.stack((state[...,2],state[...,3],laplacian(state[...,0],box/n)-fp,
                      laplacian(state[...,1],box/n)-fa),axis=-1)
        eps=1e-6;plus=compiled(state+eps*rhs,mode,label+'-plus');minus=compiled(state-eps*rhs,mode,label+'-minus')
        for i,kind in enumerate(('fixed','moving')):
            # The 21st C++ moving value is the FV enclosure radius, absent in
            # the independent energy/current reference. Do not count it as checked.
            count=len(reference_arrays[i]);scale=values[i][0]
            agreement=float(np.max(np.abs(values[i][:count]-reference_arrays[i]))/scale)
            v=values[i];expected=np.array([-v[3]-v[4],-v[3]-v[5],-v[4]+v[5]])
            derivative=(plus[i][:3]-minus[i][:3])/(2*eps)
            continuity=float(np.max(np.abs(derivative-expected))/scale)
            results.append(dict(model=label,monitor=kind,compared_values=count,
                reference_scaled_error=agreement,continuity_scaled_error=continuity,
                nonzero_exchange=bool(abs(v[5])>1e-5),
                passed=bool(agreement<gates['reference_scaled_error'] and continuity<gates['continuity_scaled_error'] and abs(v[5])>1e-5)))
    sources=['formation_monitor_probe.cpp','pocket_potential.h','formation_potential.h','pocket_monitor.h',
             'tracking_monitor.h','accurate_sum.h','monitor_reference.py','tracking_monitor_reference.py',
             'paper_reference.py','model.py','pocket_geometry.py','validate_formation_monitors.py']
    report=dict(kind=__doc__,n=n,box=box,eps=eps,parameters=dict(lambda_=8.,delta_v=.01,mass2=1.,formation=vars(Parameters())),
        normalization='absolute initial monitor total energy, for all output and continuity comparisons',
        gates=gates,results=results,passed=all(x['passed'] for x in results),
        provenance=dict(build_command=command,source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sources},
                        executable_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()))
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(passed=report['passed'],results=results),indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()

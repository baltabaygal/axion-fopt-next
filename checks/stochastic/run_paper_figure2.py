"""Evolve prepared Fig.2 physical catalogue/maps; see PAPER_FIGURE2.md."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from run_field import ROOT,read_state
from axion_fopt_next.interchange import write_arrival_map
from plot_energy import baseline,energies


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'outputs/stochastic/paper-figure2')
    parser.add_argument('--output',type=Path,default=ROOT/'checks/stochastic/results/paper_figure2')
    args=parser.parse_args();root=args.source.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    catalogue=json.loads((root/'catalogue.json').read_text())
    if catalogue['first_birth']<=10:raise AssertionError('Initial time must precede every bubble')
    c=.5;origin=5.;duration=26.625;box=4*np.pi;exe=ROOT/'build/cosmolattice/axion_fopt_evolve'
    runs=[];gates=[];histories={}
    def gate(name,value,limit):gates.append(dict(name=name,value=float(value),limit=limit,passed=bool(np.isfinite(value) and value<=limit)))
    for n,dt in ((128,.005),(64,.005),(64,.0025)):
        name=f'N{n}_dt{dt:g}';arrival=c*np.load(root/f'arrival{n}.npy');mapfile=root/f'N{n}.afopt'
        if not mapfile.exists():write_arrival_map(mapfile,arrival,box)
        directory=root/name
        params=dict(N=n,lSide=box,dt=dt,tMax=duration,evolver='PV2',expansion='true',fixedBackground='true',H0=.2,
                    omegaEoS=1/3,tau_origin=origin,omega_star=1.,f_axion=1.,theta_initial=2.,theta_prime_initial=0.,
                    mass_false=0.,mass_true=1.,potential='cosine',wall_width=0.,profile_wall_speed=.15,
                    arrival_map=str(mapfile),result_dir=str(directory),snapshot_format='binary',snapshot_at=20.,
                    output_every=round(.05/dt),energy_diagnostics='true',profile='true')
        (root/f'{name}.input.json').write_text(json.dumps(params,indent=2)+'\n')
        print(f'EVOLVE {name}: {round(duration/dt)} steps',flush=True)
        begin=time.perf_counter()
        with (root/f'{name}.log').open('w') as log:
            subprocess.run([str(exe),*[f'{k}={v}' for k,v in params.items()]],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
        seconds=time.perf_counter()-begin
        history=np.genfromtxt(directory/'history.csv',delimiter=',',names=True)
        if not all(np.all(np.isfinite(history[k])) for k in history.dtype.names):raise AssertionError('Nonfinite energy history')
        physical_tau=history['tau']/c;physical_a=c*history['a'];t=.025*physical_tau**2
        gate(name+'/scale',np.max(np.abs(physical_a-.05*physical_tau)),1e-11)
        for file in ('initial.bin','checkpoint.bin','final.bin'):
            meta,state=read_state(directory/file);program_tau=origin+meta['time']
            np.testing.assert_array_equal(state[...,2],arrival)
            np.testing.assert_array_equal(state[...,3],(arrival<=program_tau).astype(float))
            idx=np.argmin(np.abs(history['u']-meta['time']))
            gate(name+'/'+file+'/time',abs(history['u'][idx]-meta['time']),1e-11)
            k,g,v=energies(state,meta['scale'],box)
            for key,value in (('kinetic_mean',k),('gradient_mean',g),('potential_mean',v),('total_mean',k+g+v)):
                gate(name+'/'+file+'/'+key,abs(history[key][idx]-value),1e-10)
            if file=='initial.bin':
                np.testing.assert_array_equal(state[...,0],2.);np.testing.assert_array_equal(state[...,1],0.);np.testing.assert_array_equal(state[...,3],0.)
        run=dict(name=name,n=n,dt=dt,physical_dtau=dt/c,seconds=seconds,steps=round(duration/dt),
                 final_energy=float(history['total_mean'][-1]),peak_energy=float(history['total_mean'].max()),
                 peak_time=float(t[np.argmax(history['total_mean'])]),timings=json.loads((directory/'timings.json').read_text()))
        runs.append(run);histories[n,dt]=(t,history)
        np.savetxt(out/f'{name}_energy.csv',np.column_stack([t,physical_tau,physical_a,history['kinetic_mean'],history['gradient_mean'],
                       history['potential_mean'],history['total_mean'],history['inside_mean']]),delimiter=',',
                   header='m_phi_t,tau_physical,a_physical,kinetic_over_f2m2,gradient_over_f2m2,potential_over_f2m2,total_over_f2m2,true_vacuum_energy_over_f2m2',comments='')
        print(f'DONE {name}: {seconds:.2f}s; energy={run["final_energy"]:.6f}',flush=True)
    fine=runs[2]['final_energy'];gate('N64/timestep_final_energy',abs(runs[1]['final_energy']/fine-1),.01)
    time_change=float(np.max(np.abs(histories[64,.005][1]['total_mean']-histories[64,.0025][1]['total_mean'])))
    spatial=abs(runs[0]['final_energy']/runs[1]['final_energy']-1)
    t,history=histories[128,.005];nofopt=baseline(t,amplitude=2.)
    figdata=np.column_stack([t,history['kinetic_mean'],history['gradient_mean'],history['potential_mean'],history['total_mean'],nofopt])
    np.savetxt(out/'paper_style_energy.csv',figdata,delimiter=',',header='m_phi_t,kinetic,gradient,potential,total_FOPT,total_no_FOPT',comments='')
    sources=[Path(__file__),Path(__file__).with_name('PAPER_FIGURE2.md'),Path(__file__).with_name('arrival_projection.py'),
             Path(__file__).with_name('spatial_engulfment.py'),Path(__file__).with_name('tracer_reference.py'),
             Path(__file__).with_name('plot_energy.py'),ROOT/'adapters/cosmolattice/evolve.cpp',ROOT/'adapters/cosmolattice/axion_fopt.h']
    report=dict(status='passed' if all(g['passed'] for g in gates) else 'failed',parameters=dict(h=.05,beta=.2,v=.15,theta_i=2.,physical_box=8*np.pi,n=128,wall_width=0.,seed=0,tracers=100000),
                coordinate_gauge=dict(c=c,program_box=box,physical_tau_initial=10.,physical_t_initial=2.5,program_tau_initial=5.,native_H0=.2),
                runs=runs,gates=gates,catalogue=catalogue,projection=json.loads((root/'projection.json').read_text()),
                spatial_final_energy_relative_change=float(spatial),timestep_max_absolute_energy_change=time_change,
                no_FOPT_final_energy=float(nofopt[-1]),executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                limits='Fig.2 physical parameters matched. New seed realization and solver; same published-model tracers with fixed dtau=.01; not bitwise replay of original paper. Finite tracer population and spatial convergence remain diagnostics.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    for name in report['source_sha256']:
        dst=root/'sources'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dst)
    plot(out,histories,nofopt)
    print(json.dumps(dict(status=report['status'],spatial_final_energy_change=spatial,failed=[g for g in gates if not g['passed']]),indent=2))


def plot(out,histories,nofopt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,5.6),layout='constrained')
    t,h=histories[128,.005]
    # Before the numerical start the sharp-mass, pre-birth field is exactly frozen.
    for key,label,color in (('kinetic_mean','Kinetic','#ce3d3d'),('gradient_mean','Gradient','#dba521'),
                            ('potential_mean','Potential','#326ab5'),('total_mean','FOPT total','#171717')):
        ax.plot(np.r_[0,t],np.r_[0,h[key]],label=label,color=color,lw=2 if key=='total_mean' else 1.5)
    ax.plot(t,nofopt,'--',color='#777777',label='No FOPT',lw=1.5)
    ax.set(xlabel=r'$m_\phi t$',ylabel=r'$\bar\rho_\phi/(f_\phi^2 m_\phi^2)$',xlim=(0,100),ylim=(0,max(.1,float(h['total_mean'].max())*1.15)))
    ax.set_title(r'$\beta/H_*=4,\ H_*/m_\phi=0.05,\ v_w=0.15,\ \theta_i=2$'+'\nL=8π/mφ · N=128 · sharp wall · seed 0')
    ax.legend(fontsize=10);ax.grid(alpha=.15)
    for ext in ('png','pdf'):fig.savefig(out/f'energy_evolution.{ext}',dpi=200)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    for (n,dt),(time,h) in histories.items():
        ax.plot(time,h['total_mean'],label=f'N={n}, physical dτ={dt/.5:g}')
    ax.set(xlabel=r'$m_\phi t$',ylabel=r'$\bar\rho_\phi/(f_\phi^2 m_\phi^2)$',title='Fixed catalogue: grid and field-timestep comparison');ax.legend();ax.grid(alpha=.15)
    for ext in ('png','pdf'):fig.savefig(out/f'energy_refinement.{ext}',dpi=200)
    plt.close(fig)


if __name__=='__main__':main()

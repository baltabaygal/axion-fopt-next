"""First stochastic geometry plus cosine axion evolution; FIELD_EXPERIMENT.md."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import time
import numpy as np
from scipy.integrate import solve_ivp
from run_survival import ROOT,curve
from analytic import Case,survival
from tracer_reference import DEFAULT_SOURCE,load_paper,simulate
from spatial_engulfment import accelerated_paper
from axion_fopt_next.interchange import write_arrival_map

spec=importlib.util.spec_from_file_location('bubble_state_reader',ROOT/'checks/single_bubble/reference.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
read_state=reader.read_state


def mass(tau,arrival):return .5*(1+np.tanh((tau-arrival)*.5/.15))


def independent(arrival):
    n=arrival.shape[0];count=arrival.size;dx=3/n
    def rhs(u,state):
        a=1+u;theta=state[:count].reshape(arrival.shape);pi=state[count:]
        lap=sum(np.roll(theta,1,axis)+np.roll(theta,-1,axis)-2*theta for axis in range(3))/dx**2
        return np.r_[pi/a**2,(a**2*lap-a**4*mass(a,arrival)*np.sin(theta)).ravel()]
    sol=solve_ivp(rhs,(0,1.5),np.r_[np.ones(count),np.zeros(count)],method='DOP853',
                  rtol=2e-10,atol=2e-12,max_step=.01,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'outputs/stochastic/field-reference')
    parser.add_argument('--executable',type=Path,default=ROOT/'build/cosmolattice/axion_fopt_evolve')
    parser.add_argument('--paper-source',type=Path,default=DEFAULT_SOURCE)
    args=parser.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    executable=args.executable.resolve()
    shutil.copy2(args.paper_source,root/'paper_axion_sim.py.txt')
    paper,provenance=load_paper(args.paper_source)
    import numba
    numba.set_num_threads(min(8,numba.config.NUMBA_NUM_THREADS))
    case=Case('field',3.,.5,2.5);times=np.linspace(0,2.5,26)
    simulate(accelerated_paper(paper),case,92843,512,.05,times)
    begin=time.perf_counter();catalogue,internal=simulate(accelerated_paper(paper),case,0,100000,.05,times)
    catalogue_seconds=time.perf_counter()-begin
    baseline,baseline_F=simulate(paper,case,0,100000,.05,times)
    np.testing.assert_array_equal(catalogue.times,baseline.times);np.testing.assert_array_equal(catalogue.positions,baseline.positions)
    np.testing.assert_array_equal(internal,baseline_F)
    probes=np.random.default_rng(92847).uniform(0,3,(16384,3))
    geometric=curve(catalogue.arrival_times(probes),times);analytic=survival(case,times)
    np.savez_compressed(root/'catalogue.npz',times=catalogue.times,positions=catalogue.positions,internal_F=internal,
                        tau=times,geometric_F=geometric,analytic_F=analytic)
    gates=[];records=[];maps={};states={}
    def gate(name,value,limit):
        gates.append(dict(name=name,value=float(value),limit=float(limit),passed=bool(np.isfinite(value) and value<=limit)))
    for n in (8,32,64):
        print(f'PROJECT N={n}, bubbles={len(catalogue.times)}',flush=True)
        begin=time.perf_counter();maps[n]=catalogue.lattice_arrival_times(n)
        seconds=time.perf_counter()-begin
        write_arrival_map(root/f'N{n}.afopt',maps[n],3.)
        np.savetxt(root/f'N{n}_geometry.csv',np.column_stack([times,curve(maps[n].ravel(),times)]),delimiter=',',header='tau,geometric_F',comments='')
        print(f'PROJECT DONE N={n}: {seconds:.3f}s',flush=True)
        if n==8:
            begin=time.perf_counter();reference=independent(maps[n]);reference_seconds=time.perf_counter()-begin
            np.savez_compressed(root/'N8_reference.npz',u=np.array([0.,.75,1.5]),state=reference.sol([0.,.75,1.5]))
        for dt in ((.005,.0025) if n in (8,64) else (.005,)):
            name=f'N{n}_dt{dt:g}';out=root/name
            params=dict(N=n,lSide=3.,dt=dt,tMax=1.5,evolver='PV2',expansion='true',fixedBackground='true',H0=1.,
                        omegaEoS=1/3,tau_origin=1.,omega_star=1.,f_axion=1.,theta_initial=1.,theta_prime_initial=0.,
                        mass_false=0.,mass_true=1.,potential='cosine',wall_width=.15,profile_wall_speed=.5,
                        arrival_map=str(root/f'N{n}.afopt'),result_dir=str(out),snapshot_format='binary',snapshot_at=.75,
                        output_every=round(.025/dt),profile='true')
            (root/f'{name}.input.json').write_text(json.dumps(params,indent=2)+'\n')
            print(f'EVOLVE {name}: {round(1.5/dt)} steps',flush=True)
            begin=time.perf_counter()
            with (root/f'{name}.log').open('w') as log:
                subprocess.run([str(executable),*[f'{k}={v}' for k,v in params.items()]],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1200)
            elapsed=time.perf_counter()-begin
            history=np.genfromtxt(out/'history.csv',delimiter=',',names=True)
            if not all(np.all(np.isfinite(history[k])) for k in history.dtype.names):raise AssertionError('Nonfinite history')
            gate(name+'/scale_history',np.max(np.abs(history['a']-(1+history['u']))),1e-12)
            gate(name+'/tau_history',np.max(np.abs(history['tau']-(1+history['u']))),1e-12)
            reference_theta=0.;reference_pi=0.;snapstats=[]
            for file in ('initial.bin','checkpoint.bin','final.bin'):
                meta,state=read_state(out/file);u=meta['time'];a=1+u
                np.testing.assert_array_equal(state[...,2],maps[n])
                gate(name+'/'+file+'/mass',np.max(np.abs(state[...,3]-mass(a,maps[n]))),3e-14)
                gate(name+'/'+file+'/scale',abs(meta['scale']-a),1e-12)
                idx=np.argmin(np.abs(history['u']-u));gate(name+'/'+file+'/time',abs(history['u'][idx]-u),1e-12)
                field=state[...,0];pi=state[...,1];potential=np.mean(state[...,3]*(1-np.cos(field)))
                for col,value in (('theta_mean',field.mean()),('theta2_mean',np.mean(field**2)),('pi_mean',pi.mean()),('potential_mean',potential)):
                    gate(name+'/'+file+'/'+col,abs(history[col][idx]-value),1e-12)
                grad=.5*sum(np.mean(((np.roll(field,-1,axis)-field)/(3/n))**2) for axis in range(3))/a**2
                snapstats.append(dict(u=u,tau=a,t=a**2/2,theta_mean=float(field.mean()),theta_std=float(field.std()),
                                      potential=float(potential),kinetic=float(np.mean(pi**2)/(2*a**6)),gradient=float(grad),
                                      geometric_F=float(np.mean(maps[n]>a)),mean_mass_squared=float(state[...,3].mean())))
                if n==8:
                    ref=reference.sol(u);reference_theta=max(reference_theta,float(np.max(np.abs(field.ravel()-ref[:n**3]))))
                    reference_pi=max(reference_pi,float(np.max(np.abs(pi.ravel()-ref[n**3:]))))
                if file=='initial.bin':
                    np.testing.assert_array_equal(field,np.ones_like(field));np.testing.assert_array_equal(pi,np.zeros_like(pi))
                if file=='final.bin':states[n,dt]=field.copy()
            if n==8:
                gate(name+'/reference_theta',reference_theta,2e-4);gate(name+'/reference_pi',reference_pi,1e-3)
            record=dict(name=name,n=n,dt=dt,steps=round(1.5/dt),arrival_projection_seconds=seconds,wall_seconds=elapsed,
                        native_timings=json.loads((out/'timings.json').read_text()),snapshots=snapstats,
                        reference_theta_error=reference_theta if n==8 else None,reference_pi_error=reference_pi if n==8 else None)
            records.append(record);print(f'DONE {name}: {elapsed:.3f}s; theta_std={snapstats[-1]["theta_std"]:.6f}',flush=True)
    small=[r for r in records if r['n']==8]
    gate('N8/reference_refinement',3*small[1]['reference_theta_error']/small[0]['reference_theta_error'],1.)
    gate('N64/timestep_theta_rms',np.sqrt(np.mean((states[64,.005]-states[64,.0025])**2)),.001)
    averaged=states[64,.005].reshape(32,2,32,2,32,2).mean(axis=(1,3,5))
    gate('N32_N64/spatial_theta_rms',np.sqrt(np.mean((states[32,.005]-averaged)**2)),.02)
    sources=[Path(__file__),Path(__file__).with_name('FIELD_EXPERIMENT.md'),Path(__file__).with_name('spatial_engulfment.py'),
             Path(__file__).with_name('tracer_reference.py'),Path(__file__).with_name('analytic.py'),ROOT/'src/axion_fopt_next/core.py',
             ROOT/'src/axion_fopt_next/interchange.py',ROOT/'checks/single_bubble/reference.py',
             *[ROOT/'adapters/cosmolattice'/p for p in ('axion_fopt.h','potential_functions.h','evolve.cpp','prepare_overlay.py','CMakeLists.txt')]]
    report=dict(status='passed' if all(g['passed'] for g in gates) else 'failed',runs=records,gates=gates,
                catalogue=dict(bubbles=len(catalogue.times),seed=0,tracers=100000,seconds=catalogue_seconds,baseline_exact_match=True,
                               tau=times.tolist(),internal_F=internal.tolist(),geometric_F=geometric.tolist(),analytic_F=analytic.tolist()),
                reference_seconds=reference_seconds,paper_source=provenance,executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                build=json.loads((executable.parent/'build-info.json').read_text()),numpy=np.__version__,
                source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                limits='Single seeded prescribed-bubble cosine field run. Finite initial tau=1. No backreaction or abundance claim. Smooth first-arrival profile differs from sharp geometric F.')
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    for source in report['source_sha256']:
        dst=root/'sources'/source;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/source,dst)
    plot(root,report)
    print(json.dumps(dict(status=report['status'],failed=[g for g in gates if not g['passed']]),indent=2))


def plot(root,report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for r in report['runs']:
        if r['n']==8:continue
        h=np.genfromtxt(root/r['name']/'history.csv',delimiter=',',names=True);t=h['tau']**2/2
        label=f'N={r["n"]}, du={r["dt"]:g}'
        axes[0,0].plot(t,h['theta_mean'],label=label)
        axes[0,1].plot(t,np.sqrt(np.maximum(0,h['theta2_mean']-h['theta_mean']**2)),label=label)
    c=report['catalogue'];t=np.array(c['tau'])**2/2
    axes[1,0].plot(t,c['analytic_F'],color='black',label='Analytic ensemble F')
    axes[1,0].plot(t,c['geometric_F'],ls='--',label='Catalogue geometry, seed 0')
    axes[1,0].plot(t,c['internal_F'],ls=':',label='Internal tracer fraction')
    axes[1,0].axvspan(0,.5,color='grey',alpha=.1,label='Before field initial data')
    _,state=read_state(root/'N64_dt0.0025'/'final.bin')
    im=axes[1,1].imshow(state[:,:,32,0].T,origin='lower',extent=(0,3,0,3),cmap='viridis')
    fig.colorbar(im,ax=axes[1,1],label='theta');axes[1,1].set_xlabel('x (comoving)');axes[1,1].set_ylabel('y (comoving)')
    axes[1,1].set_title('Final field slice, z=1.5234, tau=2.5')
    for ax,title,y in ((axes[0,0],'Axion mean','Mean theta'),(axes[0,1],'Generated inhomogeneity','Standard deviation of theta'),
                       (axes[1,0],'False-vacuum geometry','False-vacuum fraction F')):
        ax.set_title(title);ax.set_xlabel('Cosmic time t=tau²/2');ax.set_ylabel(y);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Stochastic bubbles driving a cosine axion field\nH*=beta=mt=1 · mf=0 · L=3 · v=0.5 · 100,000 tracers · seed 0')
    for ext in ('png','pdf'):fig.savefig(root/f'field_evolution.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':main()

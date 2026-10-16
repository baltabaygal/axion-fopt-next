"""Paper-normalized physical energy from actual stochastic field snapshots."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from scipy.integrate import solve_ivp,cumulative_trapezoid
from run_field import ROOT,read_state


def energies(state,scale,box):
    """pi=a^2 theta_tau, hence theta_t=pi/a^3 (omega_star=f=1)."""
    theta=state[...,0];pi=state[...,1];m2=state[...,3];dx=box/theta.shape[0]
    kinetic=.5*np.mean(pi**2)/scale**6
    gradient=.5*sum(np.mean(((np.roll(theta,-1,axis)-theta)/dx)**2) for axis in range(3))/scale**2
    potential=np.mean(m2*(1-np.cos(theta)))
    return float(kinetic),float(gradient),float(potential)


def baseline(times,amplitude=1.,mass=1.):
    """No FOPT: regular early-time cosine misalignment in radiation domination."""
    tiny=1e-5
    initial=[amplitude-mass**2*np.sin(amplitude)*tiny**2/5,-2*mass**2*np.sin(amplitude)*tiny/5]
    def rhs(t,state):return [state[1],-3*state[1]/(2*t)-mass**2*np.sin(state[0])]
    solution=solve_ivp(rhs,(tiny,float(np.max(times))),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.1,dense_output=True)
    if not solution.success:raise RuntimeError(solution.message)
    theta,velocity=solution.sol(np.maximum(times,tiny))
    return (.5*velocity**2+mass**2*(1-np.cos(theta)))/mass**2


def process(source):
    cfg=json.loads((source/'animation.json').read_text());h=cfg['h'];box=cfg['box'];width=cfg['width']
    history=np.genfromtxt(source/'evolution/history.csv',delimiter=',',names=True)
    rows=[];max_potential=0.;max_parseval=0.
    files=sorted((source/'evolution').glob('frame_*.bin'))
    hashes=[]
    for index,path in enumerate(files):
        meta,state=read_state(path);tau=cfg['origin']+meta['time'];a=meta['scale'];t=h*tau**2/2
        np.testing.assert_allclose(a,h*tau,atol=1e-11,rtol=0)
        k,g,v=energies(state,a,box);total=k+g+v
        record=np.argmin(np.abs(history['tau']-tau))
        np.testing.assert_allclose(history['tau'][record],tau,atol=1e-11,rtol=0)
        diff=abs(v-history['potential_mean'][record]);max_potential=max(max_potential,diff)
        if diff>1e-12:raise AssertionError('Potential differs from native history')
        z=(tau-state[...,2])*.5/width
        dm2=.25/width*(1-np.tanh(z)**2)
        source_term=float(np.mean(dm2*(1-np.cos(state[...,0]))))
        rows.append([tau,t,a,k,g,v,total,source_term,6*h/a*k+2*h/a*g,float(np.mean(state[...,2]>tau))])
        if index%25==0 or index==len(files)-1:
            n=state.shape[0];wave=2*np.pi*np.fft.fftfreq(n,d=box/n)
            dispersion=(2*np.sin(wave*(box/n)/2)/(box/n))**2
            squared=dispersion[:,None,None]+dispersion[None,:,None]+dispersion[None,None,:]
            fourier=np.fft.fftn(state[...,0])/n**3
            fft_gradient=float(.5*np.sum(squared*np.abs(fourier)**2)/a**2)
            relative=abs(fft_gradient-g)/max(g,1e-30);max_parseval=max(max_parseval,relative)
            if relative>1e-12:raise AssertionError('Gradient Parseval failed')
        hashes.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    rows=np.array(rows)
    if np.any(~np.isfinite(rows)) or np.any(rows[:,3:7]<0):raise AssertionError('Invalid physical energy')
    # Radiation and explicit mass work: rho_tau+6(a_tau/a)K+2(a_tau/a)G=S.
    # Cadence dtau=.2 is coarse for this balance integral, so report it as a
    # diagnostic rather than imposing closed-system conservation.
    balance=rows[:,6]-rows[0,6]-cumulative_trapezoid(rows[:,7]-rows[:,8],rows[:,0],initial=0)
    return cfg,rows,dict(max_native_potential_difference=max_potential,max_gradient_parseval_relative_error=max_parseval,
                         balance_max_abs=float(np.max(np.abs(balance))),balance_final=float(balance[-1]),
                         balance_relative_to_peak=float(np.max(np.abs(balance))/rows[:,6].max()),snapshot_sha256=hashes)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'outputs/stochastic/representative-H0.05')
    parser.add_argument('--output',type=Path,default=ROOT/'checks/stochastic/results/energy_H0.05')
    args=parser.parse_args();source=args.source.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    cfg,rows,validation=process(source)
    nofopt=baseline(rows[:,1]);noearly=baseline(np.linspace(0,rows[-1,1],2001))
    np.savetxt(out/'energy.csv',np.column_stack([rows,nofopt,rows[:,2]**3*rows[:,6]]),delimiter=',',
               header='tau,m_phi_t,a,kinetic_over_f2m2,gradient_over_f2m2,potential_over_f2m2,total_over_f2m2,mass_work_per_tau,expansion_loss_per_tau,geometric_F,no_FOPT_total_over_f2m2,a3_total_over_f2m2',comments='')
    # A coarser-grid energy diagnostic, keeping the same catalogue/time window.
    coarse_history=np.genfromtxt(source/'N32_diagnostic/history.csv',delimiter=',',names=True)
    meta,coarse=read_state(source/'N32_diagnostic/final.bin');kc,gc,vc=energies(coarse,meta['scale'],cfg['box'])
    coarse_total=kc+gc+vc
    final=rows[-1,6]
    summary=dict(parameters=cfg,normalization='rho_bar/(f_phi^2*m_phi^2), x=m_phi*t (cosmic time); m_phi=f_phi=omega_star=1',
                 formulas=dict(kinetic='mean(pi^2)/(2*a^6)',gradient='mean(sum forward-periodic differences squared)/(2*a^2)',potential='mean(m2*(1-cos(theta)))'),
                 validation=validation,peak=dict(t=float(rows[np.argmax(rows[:,6]),1]),total=float(rows[:,6].max())),
                 final=dict(t=float(rows[-1,1]),kinetic=float(rows[-1,3]),gradient=float(rows[-1,4]),potential=float(rows[-1,5]),total=float(final),
                            no_FOPT_total=float(nofopt[-1]),a3_total=float(rows[-1,2]**3*final)),
                 spatial_diagnostic=dict(N32_final_total=coarse_total,N64_final_total=float(final),relative_difference=float(abs(coarse_total/final-1)),
                                         warning='Detailed field patterns are not spatially converged; energy convergence also requires further refinement.'),
                 no_FOPT='Independent homogeneous cosine ODE in RD with regular early theta=1. This is not an identical finite-t=10 initial condition.',
                 paper='https://arxiv.org/pdf/2607.01333, Eq.12 and Fig.2; publication parameters differ from this representative run.',
                 source_directory=str(source),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    shutil.copy2(source/'validation.json',out/'field_validation.json')
    plot(out,cfg,rows,nofopt,noearly)
    print(json.dumps({k:summary[k] for k in ('peak','final','spatial_diagnostic')},indent=2))


def plot(out,cfg,rows,nofopt,noearly):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9,5.6),layout='constrained')
    labels=('Kinetic','Gradient','Potential','Total (FOPT)')
    colors=('#cf3e3e','#dda624','#3267b2','#161616')
    for column,label,color in zip((3,4,5,6),labels,colors):ax.plot(rows[:,1],rows[:,column],color=color,lw=1.5 if column<6 else 2,label=label)
    ax.plot(rows[:,1],nofopt,'--',color='#777777',lw=1.5,label='No FOPT (homogeneous, early-time IC)')
    ax.set(xlabel=r'$m_\phi t$',ylabel=r'$\bar\rho_\phi/(f_\phi^2 m_\phi^2)$',xlim=(0,rows[-1,1]),ylim=(0,rows[:,6].max()*1.18))
    ax.axvspan(0,rows[0,1],color='#dbe1e6',alpha=.35,zorder=-2,label='Before simulated field initial data')
    ax.legend(loc='upper right',fontsize=9);ax.grid(alpha=.15)
    ax.set_title(r'$H_*/m_\phi=0.05,\quad \beta/H_*=1,\quad v_w=0.5,\quad\theta_i=1$'+'\nCosine potential · N=64 · L=60 · representative smooth-wall run')
    for ext in ('png','pdf'):fig.savefig(out/f'energy_evolution.{ext}',dpi=200)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for column,label,color in zip((3,4,5,6),labels,colors):axes[0].plot(rows[:,1],np.maximum(rows[:,column],1e-20),color=color,label=label)
    axes[0].plot(rows[:,1],nofopt,'--',color='#777777',label='No FOPT')
    axes[0].set(yscale='log',xlabel=r'$m_\phi t$',ylabel=r'$\bar\rho_\phi/(f_\phi^2 m_\phi^2)$',ylim=(1e-6,rows[:,6].max()*2));axes[0].legend(fontsize=8)
    axes[0].set_title('Physical energy density')
    axes[1].plot(rows[:,1],rows[:,2]**3*rows[:,6],color='#161616',label='FOPT')
    axes[1].plot(rows[:,1],rows[:,2]**3*nofopt,'--',color='#777777',label='No FOPT')
    axes[1].set(xlabel=r'$m_\phi t$',ylabel=r'$a^3\bar\rho_\phi/(f_\phi^2 m_\phi^2)$');axes[1].legend(fontsize=8)
    axes[1].set_title('Comoving energy diagnostic')
    for ax in axes:ax.grid(alpha=.15)
    for ext in ('png','pdf'):fig.savefig(out/f'energy_diagnostics.{ext}',dpi=200)
    plt.close(fig)


if __name__=='__main__':main()

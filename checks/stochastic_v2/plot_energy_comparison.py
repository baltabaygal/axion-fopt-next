"""Energy budgets and explicitly unmatched comparison to legacy single bubble.

Read-only analysis of archived CSVs; never starts or changes a simulation.
"""
import argparse
import json
from pathlib import Path
import textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_backend import sha

ROOT=Path(__file__).resolve().parents[2]
COLORS=dict(kinetic='#F0986E',gradient='#B8A037',potential='#5477C4',total='#1F2430')
STYLES=dict(kinetic='-',gradient='--',potential='-.',total='-')

def read(path,required):
    rows=np.genfromtxt(path,delimiter=',',names=True)
    if not set(required).issubset(rows.dtype.names):raise ValueError('Missing columns: '+str(path))
    if not all(np.all(np.isfinite(rows[c])) for c in rows.dtype.names):raise ValueError('Nonfinite CSV')
    return rows

def line(ax,x,y,label,color,style='-'):
    ax.plot(x,y,label=label,color=color,ls=style,lw=1.7)

def add_chart_header(fig,ax,title,subtitle):
    # Adapted shared-header skeleton; Matplotlib fallback because Seaborn is
    # unavailable in both local Python runtimes. Scientific PNG/SVG delivery.
    fig.subplots_adjust(top=.80 if fig.get_figheight()>7 else .75)
    fig.text(.095,.98,title,ha='left',va='top',fontsize=15,color='#1F2430')
    fig.text(.095,.927,textwrap.fill(subtitle,125),ha='left',va='top',fontsize=10,color='#6F768A')

def finish(fig,axes,out,name,title,subtitle):
    for ax in np.asarray(axes).ravel():
        ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.18,lw=.6)
        ax.legend(fontsize=8,frameon=False)
    add_chart_header(fig,np.asarray(axes).ravel()[0],title,subtitle)
    for ext in ('png','svg'):fig.savefig(out/f'{name}.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)

def components(ax,x,data,denominator):
    for key,values in data.items():line(ax,x,values/denominator,key.capitalize(),COLORS[key],STYLES[key])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--source',type=Path,default=ROOT/'outputs/stochastic_v2/single-bubble-long-01')
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    contract=dict(question='How do energy components evolve, and what can legacy/new curves actually compare?',
        takeaway='Axion energy is supplied by S; legacy growth is externally driven. Different models/boxes prevent attributing curve differences to S alone.',
        family='Trend',variant='component line small multiples plus total-energy contextual overlay',
        surface='static scientific PNG and SVG in the established simulation results folder',
        temporal_points='226 main saved times over 0..18, dense native legacy histories; no resampling except explicit common-time scalar lookup',
        palette='Three persistent component colors plus neutral total; distinct line styles; blue focal versus neutral/orange legacy curves',
        qa='Inspect final PNGs, verify component sums, normalization and energy-source balance',
        skill_limitations='Referenced user-context preflight file unavailable; Seaborn absent in local and bundled runtimes. Uses local evidence and Matplotlib template-style fallback, not full skill-template compliance.')
    (out/'chart_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    sources=[]
    def load(path,required):sources.append(path);return read(path,required)
    h=load(args.source/'main/native/history.csv',['t','E_A','E_S','total','K_A','G_A','V_A','K_S','G_S','U','Q_A'])
    half=load(args.source/'time-half/native/history.csv',['t','total'])
    short=load(ROOT/'checks/stochastic_v2/results/single-bubble-02/N96_history.csv',['t','E_A'])
    params=json.loads((args.source/'parameters.json').read_text());sources.append(args.source/'parameters.json')
    f=params['physical']['fa'];mass=params['physical']['mass'];box=48.;volume=box**3;den=volume*f*f*mass*mass
    new=dict(kinetic=h['K_A'],gradient=h['G_A'],potential=h['V_A'],total=h['E_A'])
    axion_sum_error=float(np.max(np.abs(h['K_A']+h['G_A']+h['V_A']-h['E_A'])))
    scalar_sum_error=float(np.max(np.abs(h['K_S']+h['G_S']+h['U']-h['E_S'])))
    if axion_sum_error>1e-8 or scalar_sum_error>1e-6:raise RuntimeError('Energy components do not add up')
    q=np.r_[0,np.cumsum((h['Q_A'][1:]+h['Q_A'][:-1])/2*np.diff(h['t']))]
    fig,axes=plt.subplots(2,2,figsize=(12,8.7));fig.subplots_adjust(left=.095,right=.97,bottom=.09,wspace=.28,hspace=.42)
    components(axes[0,0],h['t'],new,den)
    line(axes[0,0],short['t'],short['E_A']/den,'Original short run (N96)','#7A828F',':')
    axes[0,0].set(xlabel='Time t',ylabel=r'$\bar\rho_A/(f^2 M^2)$',title='Axion: kinetic + gradient + cosine potential',xlim=(0,18),ylim=(0,None))
    background=.06*volume
    components(axes[0,1],h['t'],dict(kinetic=h['K_S'],gradient=h['G_S'],potential=h['U']-background,total=h['E_S']-background),1.)
    axes[0,1].axhline(0,color='#7A828F',lw=.6)
    axes[0,1].set(xlabel='Time t',ylabel='Integrated energy (FV baseline subtracted)',title='S sector: subtract constant uniform-FV energy',xlim=(0,18))
    line(axes[1,0],h['t'],h['E_A']-h['E_A'][0],'Axion energy gain','#5477C4')
    line(axes[1,0],h['t'],-(h['E_S']-h['E_S'][0]),'S energy loss','#F0986E','--')
    line(axes[1,0],h['t'],q,'Integrated exchange Q_A','#1F2430',':')
    axes[1,0].set(xlabel='Time t',ylabel='Integrated energy change',title='Energy is transferred between the fields',xlim=(0,18))
    for hist,label,style in [(h,'dt = .04','-'),(half,'dt = .02','--')]:
        excess=abs(hist['total'][0]-background)
        line(axes[1,1],hist['t'],(hist['total']-hist['total'][0])/excess,label,'#5477C4',style)
    axes[1,1].set(xlabel='Time t',ylabel=r'$[E_{tot}(t)-E_{tot}(0)]/|E_{tot}(0)-E_{FV}|$',title='Closed-system numerical energy residual',xlim=(0,18))
    finish(fig,axes,out,'energy_evolution','Energy evolution of the dynamical single bubble',
        'Saved 3D CosmoLattice histories; L=48, N=128, M=.5, f=.2. Axion density is averaged over the entire box; scalar energies are volume-integrated.')
    old={}
    for v,refinement in [(.1,1),(.5,2),(.9,1)]:
        old[v]=load(ROOT/f'checks/single_bubble/results/N128_v{v}_ref{refinement}_history.csv',
            ['tau','kinetic_mean','gradient_mean','potential_mean','total_mean','source_mean'])
    legacy=old[.5]
    fig,axes=plt.subplots(1,2,figsize=(12,5.2));fig.subplots_adjust(left=.095,right=.97,bottom=.14,wspace=.28)
    for ax,x,vals,title in [(axes[0],mass*h['t'],new,'Dynamical S wall: cosine axion'),
        (axes[1],legacy['tau'],dict(kinetic=legacy['kinetic_mean'],gradient=legacy['gradient_mean'],potential=legacy['potential_mean'],total=legacy['total_mean']),'Legacy prescribed wall: harmonic axion, v_w=.5')]:
        components(ax,x,vals,den if ax is axes[0] else 1.)
        ax.set(xlabel='Mass-scaled time: M t (new) or m0 t (legacy)',ylabel='Box-mean axion density / (mass² × initial field²)',xlim=(0,9),ylim=(0,None),title=title)
    finish(fig,axes,out,'legacy_component_comparison','Axion-energy components in the new and legacy single-bubble setups',
        'Common mass-scaled time window 0..9; independent vertical scales. Different box sizes, wall profiles, initial radii and potentials: this is not a matched test.')
    fig,ax=plt.subplots(figsize=(9,5.4));fig.subplots_adjust(left=.12,right=.96,bottom=.14)
    line(ax,mass*h['t'],h['E_A']/den,'Dynamical S, cosine axion','#5477C4')
    palette={.1:'#7A828F',.5:'#F0986E',.9:'#464C55'};styles={.1:':',.5:'--',.9:'-.'}
    for v,a in old.items():line(ax,a['tau'],a['total_mean'],f'Legacy harmonic, imposed v_w={v}',palette[v],styles[v])
    ax.set(xlabel='Mass-scaled time: M t or m0 t',ylabel='Box-mean axion density / (mass² × initial field²)',yscale='log',xlim=(0,9),ylim=(1e-5,1))
    finish(fig,[ax],out,'legacy_total_comparison','Saved axion-energy histories: dynamical versus prescribed bubble walls',
        'Logarithmic density scale. New: M L=24, M R0=4; legacy: m0 L=2π, nominal R0=0. Normalization does not remove these physical differences.')
    def work_residual(a):
        work=np.r_[0,np.cumsum((a['source_mean'][1:]+a['source_mean'][:-1])/2*np.diff(a['tau']))]
        return float(np.max(np.abs(a['total_mean']-a['total_mean'][0]-work))/np.max(a['total_mean']))
    new_table=np.column_stack((h['t'],mass*h['t'],*[new[k]/den for k in ['kinetic','gradient','potential','total']],h['K_S'],h['G_S'],h['U'],h['E_S'],h['total'],h['Q_A'],q))
    np.savetxt(out/'new_energy.csv',new_table,delimiter=',',header='t,Mt,K_A_density_normalized,G_A_density_normalized,V_A_density_normalized,total_A_density_normalized,K_S_integrated,G_S_integrated,U_integrated,E_S_integrated,E_total_integrated,Q_A_integrated,cumulative_axion_exchange',comments='')
    for v,a in old.items():
        np.savetxt(out/f'legacy_v{v}_energy.csv',np.column_stack([a[n] for n in ['tau','kinetic_mean','gradient_mean','potential_mean','total_mean','source_mean']]),delimiter=',',header='m0t,K_density_normalized,G_density_normalized,V_density_normalized,total_density_normalized,external_mass_work_rate_normalized',comments='')
    audit=args.source.parent/'long-diagnostic-replay-01/report.json'
    if audit.exists():sources.append(audit)
    summary=dict(chart_contract=contract,normalization=dict(new='E_A/[L^3*f^2*M^2], theta_i=1, A_i=f',legacy='rho_phi/[m0^2*phi0^2], m0=phi0=1',new_denominator=den),
        checks=dict(axion_component_sum_max_abs=axion_sum_error,scalar_component_sum_max_abs=scalar_sum_error,
            legacy_work_balance={str(v):work_residual(a) for v,a in old.items()}),
        new_final=dict(t=float(h['t'][-1]),Mt=float(mass*h['t'][-1]),**{k:float(val[-1]/den) for k,val in new.items()},E_A=float(h['E_A'][-1]),energy_gain=float(h['E_A'][-1]-h['E_A'][0])),
        legacy_common_time={str(v):dict(m0t=min(9,float(a['tau'][-1])),total_density=float(np.interp(min(9,float(a['tau'][-1])),a['tau'],a['total_mean']))) for v,a in old.items()},
        caveat='Unmatched physics and conversion geometry; no suppression/enhancement ratio or causal attribution to adding S is inferred.',
        source_sha256={str(path.resolve()):sha(path) for path in sources},script_sha256=sha(__file__))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['checks','new_final','legacy_common_time']},indent=2))

if __name__=='__main__':main()

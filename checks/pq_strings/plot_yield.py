"""Scientific audit figure: early yield, spectral coverage, wave test and grid bound."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BLUE='#5477C4';ORANGE='#CC6F47';NEUTRAL='#464C55'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.facecolor':'white'})

def main():
    p=argparse.ArgumentParser();p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--feasibility',type=Path,required=True);p.add_argument('--wave',type=Path,required=True)
    p.add_argument('--preparation',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    audit=json.loads((args.audit/'summary.json').read_text())
    feasibility=json.loads((args.feasibility/'feasibility.json').read_text())
    wave=json.loads((args.wave/'validation.json').read_text())
    prep=json.loads((args.preparation/'summary.json').read_text())
    (args.output/'chart_contract.json').write_text(json.dumps(dict(
        surface='standalone scientific PNG/PDF',question='Can the current estimator and grids support a paper comparison?',
        panels=['early diagnostic yield line (13 snapshots each, no extrapolation)',
                'finite-box spectrum lines and unavailable IR range',
                'evolved-wave analytic normalization error lines',
                'necessary grid bound bars by mass-growth law'],
        warning='Early diagnostics are not relic yields. Bar bounds depend on stated resolution criteria.',
        palette='blue/neutral controls; orange alternate model indices',
        source_paths={k:str(getattr(args,k).resolve()) for k in ('audit','feasibility','wave','preparation')}
    ),indent=2)+'\n')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    fig.suptitle('Axion-yield estimator and feasibility audit',fontsize=16)
    fig.text(.5,.935,'Native PQ snapshots; validated isolated wave; explicit uniform-grid scale criteria',ha='center',color='#6F768A')
    for key,label,color,style in [('run_00_network_growing','Growing mass',BLUE,'-'),('run_01_network_massless','Massless control',NEUTRAL,'--')]:
        rows=audit['sources'][key]['rows']
        axes[0,0].plot([r['tau_ratio'] for r in rows],[r['continuum']['nu_diagnostic'] for r in rows],color=color,ls=style,label=label,lw=1.7)
    axes[0,0].axvline(1,color=NEUTRAL,ls=':',lw=1)
    axes[0,0].set(title='Early current diagnostic; no final yield',xlabel='tau / tau_star',ylabel='nu_ax (nonzero modes only)')
    axes[0,0].legend(frameon=False,fontsize=9)
    axes[0,0].text(.04,.93,'Dotted line: mass–H equality; late epoch ≳3 is outside view',transform=axes[0,0].transAxes,va='top',fontsize=9)
    source=args.audit/'run_00_network_growing'
    rows=audit['sources']['run_00_network_growing']['rows']
    chosen=[rows[0],rows[len(rows)//2],rows[-1]]
    for r,color,style in zip(chosen,(BLUE,ORANGE,NEUTRAL),('-','--',':')):
        data=np.genfromtxt(source/(Path(r['snapshot']).stem+'_continuum.csv'),delimiter=',',names=True)
        axes[0,1].plot(data['k_tau_star'],data['nu_per_dlnk'],color=color,ls=style,lw=1.7,label=f"tau / tau_star={r['tau_ratio']:.2f}")
    ir=rows[-1]['kmin_tau_star']
    axes[0,1].axvspan(1,ir,color='#E6E8F0',alpha=.6)
    axes[0,1].axvline(3,color=NEUTRAL,ls=':',lw=1)
    axes[0,1].set_xscale('log');axes[0,1].set_yscale('log')
    axes[0,1].set_xlim(1,420);axes[0,1].set(title='The paper’s k tau_star ≈ 3 peak is unavailable',xlabel='k tau_star',ylabel='d nu / d ln k (exact shell weights)')
    axes[0,1].legend(frameon=False,fontsize=8)
    axes[0,1].text(.04,.15,'Unavailable IR',transform=axes[0,1].transAxes,fontsize=9,color=NEUTRAL)
    for dt,color,style in [('0.005',BLUE,'-'),('0.0025',NEUTRAL,'--')]:
        history=wave['details'][dt]['history']
        axes[1,0].plot([r['tau_ratio'] for r in history],[(r['number']/r['expected_number']-1)*100 for r in history],
                       color=color,ls=style,lw=1.5,label=f'dt={dt}')
    axes[1,0].axhline(0,color=NEUTRAL,ls=':',lw=1)
    axes[1,0].set(title='Evolved wave agrees with independent linear ODE',xlabel='tau / tau_star',ylabel='Number-estimator error (%)')
    axes[1,0].legend(frameon=False,fontsize=9)
    indices=feasibility['mass_growth_indices'];ns=[r['n_tau'] for r in indices];required=[r['minimum_N'] for r in indices]
    bars=axes[1,1].bar([str(int(n)) for n in ns],required,color=[BLUE if n==8 else ORANGE for n in ns],width=.65)
    axes[1,1].set_yscale('log');axes[1,1].set_ylim(50,100000)
    axes[1,1].axhline(128,color=NEUTRAL,ls='--',lw=1,label='128 cells/side')
    for bar,n in zip(bars,required):axes[1,1].text(bar.get_x()+bar.get_width()/2,n*1.15,f'{n:,}',ha='center',fontsize=9)
    axes[1,1].set(title='Necessary grid bound for tau / tau_star = 3',xlabel='Susceptibility exponent n_tau',ylabel='Minimum N per side')
    axes[1,1].legend(frameon=False,fontsize=8,loc='upper left')
    for ax in axes.flat:ax.grid(axis='y',color='#E6E8F0',lw=.6);ax.set_axisbelow(True)
    fig.get_layout_engine().set(rect=(0,0,1,.90))
    for ext in ('png','pdf'):fig.savefig(args.output/f'yield_audit.{ext}',dpi=160)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    fig.suptitle('Core-growth preparation: separate from physical production',fontsize=14)
    for label,color,style in [('N64',BLUE,'-'),('N64_fine_dt',NEUTRAL,'--')]:
        r=prep['runs'][label];h=np.genfromtxt(args.preparation/label/'native/history.csv',delimiter=',',names=True)
        width=1/(h['a']*np.sqrt(2*h['lambda']))/(r['parameters']['lSide']/r['parameters']['N'])
        axes[0].plot(h['tau'],width,color=color,ls=style,lw=1.7,label=label)
        axes[1].plot([s['tau'] for s in r['snapshots']],[s['plaquettes'] for s in r['snapshots']],color=color,ls=style,lw=1.7)
    for ax in axes:ax.axvline(12,color=ORANGE,ls=':',lw=1.4);ax.set_xlabel('Conformal time tau');ax.grid(axis='y',color='#E6E8F0')
    axes[0].set(title='Core width grows linearly until tau_cg=12',ylabel='Characteristic core width / dx');axes[0].legend(frameon=False)
    axes[1].set(title='The larger-box network survives the handoff',ylabel='String-pierced plaquettes')
    for ext in ('png','pdf'):fig.savefig(args.output/f'preparation.{ext}',dpi=160)
    print(args.output)
if __name__=='__main__':main()

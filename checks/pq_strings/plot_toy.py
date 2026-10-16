"""Standalone scientific figures for the explicitly modified finite-box toy."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);a=p.parse_args()
    r=json.loads((a.source/'comparison.json').read_text())
    data={k:json.loads((a.source/f'{k}_toy_analysis.json').read_text()) for k in ('base','time','space','seed0','seed1','seed2')}
    control=json.loads((a.source/'control_topology.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    fig.suptitle('Capped-mass PQ toy: collapse and late angular yield',fontsize=16)
    axes[0,0].plot([x['tau'] for x in control],[x['plaquettes'] for x in control],label='Massless control',color='#555C68',ls='--')
    for key,label,col in [('base','64³, dt=0.02','#5477C4'),('space','128³, dt=0.02','#CC6F47')]:
        rows=json.loads((a.source/f'{key}_topology.json').read_text());axes[0,0].plot([x['tau'] for x in rows],[x['plaquettes'] for x in rows],label=label,color=col)
    axes[0,0].set(xlim=(40,66),xlabel='Conformal time τ',ylabel='String-pierced plaquettes (raw)',title='Resolved collapse; control retains four strings')
    # Raw counts scale with N for these axis-aligned strings; show in equivalent straight strings.
    axes[0,0].lines[0].set_ydata(np.array(axes[0,0].lines[0].get_ydata())/64)
    axes[0,0].lines[1].set_ydata(np.array(axes[0,0].lines[1].get_ydata())/64)
    axes[0,0].lines[2].set_ydata(np.array(axes[0,0].lines[2].get_ydata())/128)
    axes[0,0].set_ylabel('Plaquettes / N (straight-string equivalent)')
    axes[0,0].set_ylim(-.1,4.2)
    for key,label,style in [('base','64³, dt=0.02','-'),('time','64³, dt=0.01','--'),('space','128³, dt=0.02',':')]:
        rows=data[key]['rows'];t=[x['tau'] for x in rows];v=[np.nan if x['nu_total_lattice'] is None else x['nu_total_lattice'] for x in rows]
        axes[0,1].plot(t,v,label=label,ls=style)
    axes[0,1].set(xlim=(340,400),xlabel='Conformal time τ',ylabel='ν toy (nonzero modes + coherent mode)',ylim=(7.8,9.3),title='Late means: 6.4% change under spatial refinement')
    for i in range(3):
        rows=data[f'seed{i}']['rows'];axes[1,0].plot([x['tau'] for x in rows],[np.nan if x['nu_total_lattice'] is None else x['nu_total_lattice'] for x in rows],label=f'Seed {i}')
    axes[1,0].set(xlim=(340,440),xlabel='Conformal time τ',ylabel='ν toy',title='Random realizations: eligible frames; seed 1 extended')
    rows=data['base']['rows'];t=[x['tau'] for x in rows]
    for field,label in [('rho_total','Total PQ energy'),('rho_axion_harmonic','Harmonic angular proxy'),('rho_tilt','Tilt energy')]:
        axes[1,1].semilogy(t,[np.nan if x[field] is None else x[field] for x in rows],label=label)
    axes[1,1].set(xlim=(40,400),xlabel='Conformal time τ',ylabel='Physical energy density (f_a=1)',title='Radial energy remains separate from axion proxy')
    for ax in axes.flat:ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
    fig.savefig(a.source/'toy_results.png',dpi=180);fig.savefig(a.source/'toy_results.pdf');plt.close(fig)
if __name__=='__main__':main()

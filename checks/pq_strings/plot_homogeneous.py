"""Plot saved homogeneous comparisons without rerunning evolution."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    summary = json.loads((args.source/'summary.json').read_text())
    if not summary['passed']:
        raise ValueError('Require accepted numerical run')
    args.output.mkdir(parents=True, exist_ok=False)
    contract = dict(question='How does homogeneous PQ approach theta-only evolution as radial mass increases?',
                    family='line', variant='three-panel time series, one figure per amplitude',
                    source=str(args.source.resolve()), samples=401, time='conformal u=tau-1, 0..8',
                    normalization='phase and residual divided by fixed input amplitude; radius is r/f_a',
                    palette=['#1e3a5f','#477da8','#9dbdd3','#252525'],
                    distinction='different line styles and labelled mass ratios; neutral dashed theta-only reference',
                    exports=['PNG','PDF'], qa='inspect exported PNGs')
    (args.output/'chart_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False,
                         'axes.titleweight':'bold', 'savefig.facecolor':'white'})
    for amplitude in sorted({r['amplitude'] for r in summary['runs']}):
        runs = [r for r in summary['runs'] if r['amplitude']==amplitude]
        dt = min(r['dt'] for r in runs)
        runs = sorted((r for r in runs if r['dt']==dt), key=lambda r:r['lam'])
        fig, axes = plt.subplots(3,1,figsize=(11,9),sharex=True,layout='constrained')
        fig.suptitle(f'Homogeneous PQ field: input amplitude {amplitude:g}',fontsize=18,ha='left',x=.075)
        axes[0].set_title('Constant axion mass; radiation expansion; radial mass / axion mass = 2, 8, 32',fontsize=11,loc='left',fontweight='normal',pad=14)
        for i,(r,color,style) in enumerate(zip(runs,contract['palette'][:3],['-','--','-.'])):
            d = np.genfromtxt(args.source/r['name']/'comparison.csv',delimiter=',',names=True)
            label=f"PQ: mass ratio {r['radial_mass']:g}"
            axes[0].plot(d['u'],d['theta_pq']/amplitude,color=color,ls=style,lw=1.8,label=label)
            axes[1].plot(d['u'],(d['theta_pq']-d['theta_fixed_radius'])/amplitude,color=color,ls=style,lw=1.6)
            axes[2].plot(d['u'],d['r_pq']-1,color=color,ls=style,lw=1.4)
        axes[0].plot(d['u'],d['theta_fixed_radius']/amplitude,color='#252525',ls=':',lw=2,label='Theta-only')
        axes[0].legend(ncol=2,frameon=False,fontsize=10,loc='upper right')
        axes[0].set_ylabel(r'Phase $\theta/A$')
        axes[1].set_ylabel(r'Phase difference $(\theta_{PQ}-\theta_{only})/A$')
        axes[2].set_ylabel(r'Radial displacement $r/f_a-1$')
        axes[2].set_xlabel(r'Conformal time since start $u=\tau-1$')
        for ax in axes:
            ax.axhline(0,color='#b0b0b0',lw=.6,zorder=0)
            ax.grid(axis='y',color='#eeeeee',lw=.7)
            ax.set_xlim(0,8)
        fig.supxlabel('A = fixed input amplitude; f_a = m_a = 1. Compiled PQ solver, finest timestep 0.000025.',fontsize=10)
        for extension in ('png','pdf'):
            fig.savefig(args.output/f'homogeneous_A{amplitude:g}.{extension}',dpi=170)
        plt.close(fig)


if __name__=='__main__':
    main()

"""Plot the global conversion criterion and its numerical controls."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from expanding_model import cosmic

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('summary','output'):p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();r=json.loads(args.summary.read_text());out=args.output
    if out.exists():raise FileExistsError(out)
    fig,axes=plt.subplots(2,1,figsize=(10,8),layout='constrained',sharex=True)
    for name,label,color,style in [('main','96³ main','#245C93','-'),('space','128³ spatial control','#41865A','--'),('time_half','96³ half timestep','#CB7B30',':')]:
        case=r['results'][name]['case'];path=Path(r['source_reports'][name]['path'])/case/'native/history.csv'
        h=np.genfromtxt(path,delimiter=',',names=True)
        axes[0].plot(h['t_cosmic'],1-h['F_half'],label=label,color=color,ls=style)
        axes[1].plot(h['t_cosmic'],h['S_min'],label=label,color=color,ls=style)
        if name=='main':
            axes[0].plot(h['t_cosmic'],h['near_true_fraction'],'-.',color='gray',label='Main fraction |S−1|<0.1')
            for v in r['results'][name]['all_cell_intervals']:
                if v['duration']>=2:
                    for ax in axes:ax.axvspan(cosmic(v['start'],.025),cosmic(v['end'],.025),alpha=.12,color='#41865A')
    axes[0].axhline(1,color='black',lw=.6);axes[0].set(ylabel='Fraction of all 3D cells',ylim=(-.02,1.04),title='Phase conversion and proximity to S=1')
    axes[1].axhline(.5,color='black',lw=.8,label='Every-cell phase threshold')
    axes[1].axhline(.9,color='gray',lw=.6,ls='--');axes[1].set(xlabel='Elapsed cosmic time',ylabel='Minimum S across the full box',title='A full box requires minimum S ≥ 0.5')
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.legend(fontsize=9)
    subtitle='Shading: sustained all-cell conversion in the main run' if r['results']['main']['sustained_conversion'] else 'No sustained all-cell conversion in the main run'
    fig.suptitle('Longer stochastic evolution with prescribed radiation expansion\n'+subtitle,fontsize=13)
    fig.savefig(out,dpi=170);plt.close(fig)

if __name__=='__main__':main()

"""Static exploratory scientific plot from saved native pressure histories."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.pilot/'report.json').read_text())
    contract=dict(question='Does stronger axion field pressure reverse an initially contracting planar wall?',
        takeaway='The A=.16 pilot reverses; A=.08 slows contraction; the empty wall continues inward.',
        family='line',variant='Two aligned time-series panels, position and phase velocity',
        surface='Standalone static scientific PNG/PDF, not a dashboard',
        data=dict(series=3,points_per_series=721,time_range=[0,180],grain=.25,units='Canonical dimensionless Minkowski fields/time'),
        palette=dict(policy='hard two-root cap plus neutral',empty='#464C55',moderate='#A3BEFA',strong='#F0986E'),
        non_color='Dashed empty, dash-dot moderate, solid strong; amplitude labels',
        limits='Exploratory single-grid pilot, no statistical bands and no converged lifetime claim',
        renderer='Matplotlib scientific fallback; Seaborn unavailable and optional user-context companion missing; no template compliance claim')
    (args.output/'chart-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,'axes.labelsize':11})
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,gridspec_kw={'hspace':.28})
    fig.subplots_adjust(top=.78,bottom=.13,left=.12,right=.96)
    colors=['#464C55','#A3BEFA','#F0986E'];styles=['--','-.','-'];histories={};sources={}
    for row,color,style in zip(report['rows'],colors,styles):
        amplitude=row['amplitude'];label='empty' if amplitude==0 else f'amplitude-{amplitude}'
        path=args.pilot/label/'native/history.csv';h=np.genfromtxt(path,delimiter=',',names=True);histories[amplitude]=h
        sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        name='No axion field' if amplitude==0 else f'A = {amplitude:.2f}'
        axes[0].plot(h['t'],h['wall_x'],color=color,ls=style,lw=2.5,label=name)
        axes[1].plot(h['t'],h['wall_phase_velocity'],color=color,ls=style,lw=2.5)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.spines[['left','bottom']].set_color('#B7BCC7')
        ax.grid(axis='y',color='#E6E8F0',lw=.7)
        ax.set_xlim(0,180);ax.tick_params(colors='#464C55')
    axes[0].set_ylabel('Wall position x');axes[0].set_ylim(197,259)
    axes[1].set_ylabel('Wall phase velocity');axes[1].set_ylim(-.34,.29)
    axes[1].set_xlabel('Time t (dimensionless)');axes[1].axhline(0,color='#6F768A',lw=1.,zorder=0)
    strong=report['rows'][-1]
    axes[0].annotate('Turnaround',xy=(strong['minimum_position_time'],strong['minimum_wall_position']),
        xytext=(116,245),arrowprops={'arrowstyle':'->','color':'#464C55','lw':1.},color='#1F2430')
    fig.suptitle('Axion field pressure reverses the wall in a planar pilot',x=.12,y=.97,ha='left',fontsize=16,color='#1F2430')
    fig.text(.12,.915,'Coupled fields only · N=4096 · spatial order 4 · dt=0.0125 · initial approach speed 0.3',ha='left',fontsize=10,color='#6F768A')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper left',bbox_to_anchor=(.113,.885),ncol=3,frameon=False)
    fig.text(.12,.045,'Single-grid finite-amplitude result; refinement and boundary controls pending. Not a 3D pocket lifetime.',fontsize=9,color='#6F768A')
    fig.savefig(args.output/'wall_pressure.png',dpi=170,facecolor='white')
    fig.savefig(args.output/'wall_pressure.pdf',facecolor='white');plt.close(fig)
    sources[str(args.pilot/'report.json')]=hashlib.sha256((args.pilot/'report.json').read_bytes()).hexdigest()
    (args.output/'provenance.json').write_text(json.dumps(dict(sources=sources,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')

if __name__=='__main__':main()

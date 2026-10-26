"""Reproducible scientific overview of the completed global-stiffness pair."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    contract=dict(question='Do the grids agree, and do the completed fields satisfy the candidate screen?',
        takeaway='Histories agree closely; sector depletion and minimum wall resolution fail.',
        charts=['radius trend (zoomed)','moving-window net energy reduction trend','axion-sector energy depletion trend','final full-3D minimum wall width bars'],
        grain='321 saved time points per grid, t=0..80; one final full-3D minimum per grid',
        palette='blue/orange for N=320/384; dashed/solid redundant grid encoding',
        caveat='Static toy fields. Net window reduction is not irreversible radiation; no long-lived relic is established.',
        source=str(args.source.resolve()),sources={})
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
        'axes.spines.right':False,'axes.labelcolor':'#1F2430','text.color':'#1F2430',
        'axes.edgecolor':'#B6BDCB','grid.color':'#E6E8F0','grid.linewidth':.7,
        'figure.facecolor':'white','savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(12,8.2));fig.subplots_adjust(left=.09,right=.97,bottom=.13,top=.82,hspace=.47,wspace=.27)
    colors=['#5477C4','#CC6F47'];styles=['--','-'];widths=[]
    for name,color,style in zip(('n320','n384'),colors,styles):
        root=args.source/name;path=root/'native/history.csv'
        h=np.genfromtxt(path,delimiter=',',names=True);r=json.loads((root/'summary.json').read_text())
        assert h['u'][-1]==80 and np.isfinite(h['tracking_energy']).all()
        contract['sources'][str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        contract['sources'][str((root/'summary.json').resolve())]=hashlib.sha256((root/'summary.json').read_bytes()).hexdigest()
        label=f"N = {name[1:]}"
        axion=sum(h[k] for k in ('axion_kinetic','axion_gradient','coupling'))
        curves=(h['fv_radius'],100*(1-h['tracking_energy']/h['tracking_energy'][0]),100*(1-axion/axion[0]))
        for ax,values in zip(axes.flat,curves):
            ax.plot(h['u'],values,color=color,ls=style,lw=1.8,label=label)
        widths.append(r['full_final_wall']['min_cells'])
    titles=['A  Pocket radius (zoomed scale)','B  Net energy reduction in moving window',
            'C  Axion-sector energy depletion','D  Thinnest final 3D wall patch']
    labels=['Volume-equivalent radius','% of initial window energy','% of initial axion-sector energy','10–90% wall width / cell size']
    for ax,title,label in zip(axes.flat,titles,labels):
        ax.set_title(title,loc='left',fontsize=12,pad=12);ax.set_ylabel(label);ax.grid(axis='y',alpha=.8);ax.set_axisbelow(True)
    for ax in list(axes.flat)[:3]:ax.set_xlim(0,80);ax.set_xlabel('Time (dimensionless)');ax.set_xticks(np.arange(0,81,20))
    axes[0,0].set_ylim(18.9,20.15);axes[0,0].axhline(20,color='#6F768A',lw=1,ls=':',label='Initial radius')
    axes[0,0].legend(frameon=False,fontsize=9,loc='upper right')
    axes[0,1].set_ylim(-.05,1.4)
    axes[0,1].text(.03,.91,'Fine grid: 1.220% reduction at t = 80',transform=axes[0,1].transAxes,fontsize=10)
    axes[1,0].axhline(5,color='#464C55',lw=1,ls=':')
    axes[1,0].text(3,5.15,'Original screen: ≤5%',fontsize=9,color='#464C55');axes[1,0].set_ylim(-1,8)
    ax=axes[1,1];bars=ax.bar([0,1],widths,color=colors,width=.5,edgecolor='#464C55',linewidth=.8)
    ax.set_xticks([0,1],['N = 320','N = 384']);ax.set_ylim(0,5)
    ax.axhline(4,color='#464C55',ls=':',lw=1);ax.text(-.25,4.15,'Required: ≥4 cells',fontsize=9)
    for bar,value in zip(bars,widths):ax.text(bar.get_x()+bar.get_width()/2,value+.12,f'{value:.2f}',ha='center',fontsize=11)
    fig.text(.09,.955,'Strong confinement, but the original screens still fail',fontsize=18,weight='bold',ha='left')
    fig.text(.09,.91,'Coupled axion–dilaton fields • static background • same continuum initial state • L = 64, R₀ = 20',fontsize=11,color='#6F768A')
    fig.text(.09,.062,'Moving-window energy change includes boundary transport; axion depletion also includes exchange with the dilaton.',fontsize=10,color='#464C55')
    fig.text(.09,.035,'Two initial diameter-crossing times only. Wall resolution, time refinement, isolation and long-term stability remain open.',fontsize=10,color='#464C55')
    for extension in ('png','pdf'):fig.savefig(args.output/f'field_pair.{extension}',dpi=180)
    plt.close(fig)
    contract['plot_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.output/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')


if __name__=='__main__':main()

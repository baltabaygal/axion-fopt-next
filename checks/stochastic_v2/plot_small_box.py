"""Show actual periodic conversion, including global diagnostics and slice limits."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from model import read_state,PARAMS
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'report.json').read_text());representative='time-half'
    h=np.genfromtxt(args.source/representative/'native/history.csv',delimiter=',',names=True)
    frames=[read_state(f) for f in sorted((args.source/representative/'native').glob('frame_*.bin'))]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.2),layout='constrained')
    for name,label,style in [('main','N96, dt=.02','-'),('coarse','N64, dt=.02',':'),('time-half','N96, dt=.01','--')]:
        a=np.genfromtxt(args.source/name/'native/history.csv',delimiter=',',names=True)
        axes[0,0].plot(a['t'],1-a['F_half'],style,label=label)
    axes[0,0].axhline(1,color='gray',lw=.7);axes[0,0].set(xlabel='Time t',ylabel='Fraction of ALL 3D cells with S >= 0.5',ylim=(0,1.05),title='True-phase volume fraction');axes[0,0].legend(fontsize=8)
    axes[0,1].fill_between(h['t'],h['S_min'],h['S_max'],alpha=.18,color='#5477C4',label='Global min-to-max range')
    axes[0,1].plot(h['t'],h['S_mean'],color='#5477C4',label='Spatial mean S')
    axes[0,1].axhline(.5,color='gray',ls=':',label='Operational phase threshold');axes[0,1].axhline(1,color='gray',lw=.7)
    axes[0,1].set(xlabel='Time t',ylabel='S',title='Phase conversion does not imply S = 1 everywhere');axes[0,1].legend(fontsize=8)
    den=24**3*PARAMS['fa']**2*PARAMS['mass']**2
    for key,label,color,style in [('K_A','Kinetic','#F0986E','-'),('G_A','Gradient','#B8A037','--'),('V_A','Potential','#5477C4','-.'),('E_A','Total','#1F2430','-')]:
        axes[1,0].plot(h['t'],h[key]/den,color=color,ls=style,label=label)
    axes[1,0].set(xlabel='Time t',ylabel=r'$\bar\rho_A/(f^2 M^2)$',title='Axion energy during image collisions');axes[1,0].legend(fontsize=8)
    axes[1,1].plot(h['t'],h['near_true_fraction'],label='Fraction |S-1| < .1')
    axes[1,1].plot(h['t'],1-h['F_barrier'],'--',label='Fraction above bare barrier')
    axes[1,1].set(xlabel='Time t',ylabel='Fraction of all cells',ylim=(0,1.05),title='Proximity to the minimum is a separate test');axes[1,1].legend(fontsize=8)
    fig.suptitle('Smaller periodic box: L=24, same seeded bubble and potential\nFields and energies: N96, dt=.01; periodic-image collisions included',fontsize=14)
    fig.savefig(out/'completion.png',dpi=160);plt.close(fig)
    # Subsample only rendering cadence, never phase diagnostics used for completion.
    shown=frames[::2]
    if shown[-1][0]['t']!=frames[-1][0]['t']:shown.append(frames[-1])
    limits=[max(1.1,max(float(np.max(abs(p[...,0]))) for _,p in shown)),
        max(1.2,max(float(np.max(p[...,0]**2)) for _,p in shown)),
        max(1.1,max(float(np.max(abs(p[...,1]/PARAMS['fa']))) for _,p in shown))]
    fig,ax=plt.subplots(1,3,figsize=(12,4.4),layout='constrained');images=[]
    data=lambda y:(y[...,0],y[...,0]**2,y[...,1]/PARAMS['fa'])
    for axis,a,title,cmap,lim in zip(ax,data(shown[0][1]),['Transition field S','Axion mass squared / M²','Axion angle theta'],['viridis','magma','RdBu_r'],limits):
        lower=-lim if cmap=='RdBu_r' else (min(-.05,min(float(p[...,0].min()) for _,p in shown)) if cmap=='viridis' else 0)
        image=axis.imshow(a.T,origin='lower',extent=(-12,12,-12,12),cmap=cmap,vmin=lower,vmax=lim)
        images.append(image);axis.set(title=title,xlabel='x',ylabel='y');fig.colorbar(image,ax=axis,shrink=.7)
    heading=fig.suptitle('Periodic bubble | central slice; conversion is measured over the full 3D box')
    def update(i):
        meta,y=shown[i];j=int(np.argmin(abs(h['t']-meta['t'])))
        for image,a in zip(images,data(y)):image.set_data(a.T)
        heading.set_text(f'Periodic bubble | t={meta["t"]:.1f} | 3D true-phase fraction={1-h["F_half"][j]:.5f}')
        return images+[heading]
    ani=FuncAnimation(fig,update,frames=len(shown),interval=110,blit=False);ani.save(out/'completion.gif',writer=PillowWriter(fps=9),dpi=90);plt.close(fig)
    (out/'plot_provenance.json').write_text(json.dumps(dict(source=str(args.source.resolve()),report_sha256=sha(args.source/'report.json'),script_sha256=sha(__file__),rendered_frames=len(shown),diagnostic_cadence=.04),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()

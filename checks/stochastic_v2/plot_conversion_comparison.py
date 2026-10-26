"""Compare global conversion in L24 and L20, with refined transient checks."""
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
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--l24',type=Path,required=True)
    p.add_argument('--l20',type=Path,required=True);p.add_argument('--refinement',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    paths=dict(l24=args.l24/'time-half/native/history.csv',l20=args.l20/'main/native/history.csv',
        spatial=args.refinement/'spatial/native/history.csv',temporal=args.refinement/'temporal/native/history.csv')
    history={k:np.genfromtxt(v,delimiter=',',names=True) for k,v in paths.items()}
    base=json.loads((args.l20/'report.json').read_text());ref=json.loads((args.refinement/'report.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for name,label,color,style in [('l24','L=24, N96, dt=.01','#7A828F','--'),('l20','L=20, N80, dt=.01','#5477C4','-')]:
        h=history[name];axes[0].plot(h['t'],1-h['F_half'],color=color,ls=style,label=label)
    axes[0].axhline(1,color='gray',lw=.7);axes[0].set(xlabel='Time t',ylabel='3D volume fraction with S >= 0.5',ylim=(0,1.035),title='Smaller box reaches the threshold sooner');axes[0].legend(fontsize=8)
    for name,label,color,style in [('l20','N80, dt=.01','#5477C4','-'),('spatial','N96, dt=.01','#F0986E','--'),('temporal','N80, dt=.005','#464C55',':')]:
        h=history[name];axes[1].plot(h['t'],h['S_min'],color=color,ls=style,label=label)
    axes[1].axhline(.5,color='gray',lw=.8,label='All cells converted when min(S) >= .5')
    axes[1].set(xlabel='Time t',ylabel='Minimum S over ALL 3D cells',xlim=(27,30),ylim=(0,1.02),title='A transient all-cell interval, followed by recrossing');axes[1].legend(fontsize=8)
    fig.suptitle('Periodic single bubble: full phase conversion is temporary\nSame potential and seed; no expansion or friction',fontsize=13)
    for ext in ('png','svg'):fig.savefig(out/f'conversion_comparison.{ext}',dpi=170)
    plt.close(fig)
    # L20 original pilot with direct 3D global phase statistics in every title.
    frames=[read_state(f) for f in sorted((args.l20/'main/native').glob('frame_*.bin'))][::2]
    h=history['l20'];scalars=np.concatenate([y[...,0].ravel() for _,y in frames]);angles=np.concatenate([y[...,1].ravel()/PARAMS['fa'] for _,y in frames])
    smin=min(-.05,float(scalars.min()));smax=max(1.05,float(scalars.max()));amax=max(1.1,float(abs(angles).max()))
    fig,axes=plt.subplots(1,3,figsize=(12,4.4),layout='constrained');images=[]
    data=lambda y:(y[...,0],y[...,0]**2,y[...,1]/PARAMS['fa'])
    for ax,arr,title,cmap,lower,upper in zip(axes,data(frames[0][1]),['Transition field S','Axion mass squared / M²','Axion angle theta'],['viridis','magma','RdBu_r'],[smin,0,-amax],[smax,max(smin*smin,smax*smax),amax]):
        im=ax.imshow(arr.T,origin='lower',extent=(-10,10,-10,10),cmap=cmap,vmin=lower,vmax=upper);images.append(im)
        ax.set(title=title,xlabel='x',ylabel='y');fig.colorbar(im,ax=ax,shrink=.7)
    title=fig.suptitle('L=20 periodic box | central slices; fraction measured in all 3D cells')
    def update(i):
        meta,y=frames[i];j=int(np.argmin(abs(h['t']-meta['t'])))
        for im,a in zip(images,data(y)):im.set_data(a.T)
        title.set_text(f'L=20 | t={meta["t"]:.1f} | 3D true-phase fraction={1-h["F_half"][j]:.5f}')
        return images+[title]
    ani=FuncAnimation(fig,update,frames=len(frames),interval=120,blit=False);ani.save(out/'smaller_box.gif',writer=PillowWriter(fps=8),dpi=90);plt.close(fig)
    metrics={}
    for name,h in history.items():
        j=int(np.argmin(h['F_half']));metrics[name]=dict(best_fraction=float(1-h['F_half'][j]),best_time=float(h['t'][j]),final_fraction=float(1-h['F_half'][-1]))
    (out/'summary.json').write_text(json.dumps(dict(refinement_passed=ref['passed'],original_sustained_goal=base['goal_reached'],
        original_full_intervals=base['results']['full_intervals'],metrics=metrics,
        sources={str(v.resolve()):sha(v) for v in paths.values()},script_sha256=sha(__file__),
        note='Exact threshold fractions are grid measurements; no settled-vacuum or spatial-continuum claim'),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()

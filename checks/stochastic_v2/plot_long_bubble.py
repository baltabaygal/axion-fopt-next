"""Plot the completed long-run evidence without rerunning dynamics."""
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
    report=json.loads((args.source/'report.json').read_text());result=report['results']['main'];n=result['n'];box=result['box'];dx=box/n
    frames=[read_state(f) for f in sorted((args.source/'main/native').glob('frame_*.bin'))]
    h=np.genfromtxt(args.source/'main/native/history.csv',delimiter=',',names=True)
    extent=(-box/2,box/2,-box/2,box/2);x=(np.arange(n)+.5)*dx-box/2
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for j in (0,15,30,45):
        meta,plane=frames[j];axes[0,0].plot(x,plane[:,n//2-1:n//2+1,0].mean(axis=1),label=f't={meta["t"]:g}')
    axes[0,0].set(xlim=(0,23),xlabel='Distance from centre',ylabel='S',title='Transition-field wall');axes[0,0].legend(fontsize=9)
    for name,label,style in [('main','Axion, N128','-'),('coarse','Axion, N96',':'),('empty','Empty, N128','--')]:
        c=report['results'][name]['wall_trajectory'];t=np.array([v['t'] for v in c]);r=np.array([np.nan if v['radius'] is None else v['radius'] for v in c])
        axes[0,1].plot(t,r,style,label=label)
        if name!='coarse':axes[1,0].plot(t,np.gradient(r,t),style,label=label)
    axes[0,1].set(xlabel='Time t',ylabel='Radius at S=0.5',title='Bubble expansion');axes[0,1].legend(fontsize=9)
    axes[1,0].set(xlabel='Time t',ylabel='Measured dR/dt',title='Wall speed evolves');axes[1,0].legend(fontsize=9)
    theta_max=max(1.1,max(float(np.max(np.abs(plane[...,1]/PARAMS['fa']))) for _,plane in frames))
    im=axes[0,2].imshow((frames[-1][1][...,1]/PARAMS['fa']).T,origin='lower',extent=extent,cmap='RdBu_r',vmin=-theta_max,vmax=theta_max)
    axes[0,2].set(xlabel='x',ylabel='y',title='Axion angle at t=18');fig.colorbar(im,ax=axes[0,2],label='theta')
    q=np.r_[0,np.cumsum((h['Q_A'][1:]+h['Q_A'][:-1])/2*np.diff(h['t']))]
    axes[1,1].plot(h['t'],h['E_A']-h['E_A'][0],label='Axion energy gain')
    axes[1,1].plot(h['t'],q,'--',label='Integrated transfer from S')
    axes[1,1].plot(h['t'],-(h['E_S']-h['E_S'][0]),':',label='S energy loss')
    axes[1,1].set(xlabel='Time t',ylabel='Energy change',title='Energy exchange');axes[1,1].legend(fontsize=8)
    for name,label in [('main','With axion'),('empty','Empty')]:
        c=report['results'][name]['wall_trajectory'];axes[1,2].plot([v['t'] for v in c],[np.nan if v['width'] is None else v['width']/dx for v in c],label=label)
    axes[1,2].axhline(5,color='gray',ls=':',label='Declared minimum')
    axes[1,2].set(xlabel='Time t',ylabel='10%-90% thickness / dx',title='Saved wall-resolution diagnostic');axes[1,2].legend(fontsize=8)
    fig.suptitle('Longer single bubble: dynamical S + cosine axion | t = 0 to 18\nSame potential and seed; no expansion, friction or prescribed speed',fontsize=14)
    fig.savefig(out/'long_bubble.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    specs=[('Transition field S','viridis',-.05,1.1),('Axion mass squared / M²','magma',0,1.25),('Axion angle theta','RdBu_r',-theta_max,theta_max)]
    data=lambda plane:(plane[...,0],plane[...,0]**2,plane[...,1]/PARAMS['fa'])
    images=[]
    for axis,arr,(title,cmap,vmin,vmax) in zip(ax,data(frames[0][1]),specs):
        im=axis.imshow(arr.T,origin='lower',extent=extent,cmap=cmap,vmin=vmin,vmax=vmax);images.append(im)
        axis.set(title=title,xlabel='x',ylabel='y');fig.colorbar(im,ax=axis,shrink=.7)
    title=fig.suptitle('Dynamical bubble | t=0 | central slice of a 3D simulation')
    def update(j):
        meta,plane=frames[j]
        for image,arr in zip(images,data(plane)):image.set_data(arr.T)
        title.set_text(f'Dynamical bubble | t={meta["t"]:.1f} | central slice of a 3D simulation')
        return images+[title]
    ani=FuncAnimation(fig,update,frames=len(frames),interval=140,blit=False)
    ani.save(out/'long_bubble.gif',writer=PillowWriter(fps=7),dpi=95);plt.close(fig)
    (out/'plot_provenance.json').write_text(json.dumps(dict(source=str(args.source.resolve()),report_sha256=sha(args.source/'report.json'),renderer_sha256=sha(__file__)),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()

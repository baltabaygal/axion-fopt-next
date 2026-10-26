"""Render saved actual field evolution; no new solver runs."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from model import read_state,PARAMS

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();dest=args.output.resolve();dest.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'report.json').read_text());result=report['results']['N96']
    native=args.source/'N96/native';frames=[read_state(f) for f in sorted(native.glob('frame_*.bin'))]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    x=(np.arange(96)+.5)*.5-24
    for i in (0,len(frames)//2,len(frames)-1):
        meta,plane=frames[i];line=plane[:,47:49,0].mean(axis=1)
        ax[0,0].plot(x,line,label=f't = {meta["t"]:g}')
    ax[0,0].set(xlim=(0,16),xlabel='Distance from centre',ylabel='Transition field S',title='A resolved wall moves outward')
    ax[0,0].axhline(.5,color='gray',lw=.7,ls=':');ax[0,0].legend()
    for name,label in [('N96','With axion, N=96'),('N64','With axion, N=64'),('empty-N64','No axion, N=64')]:
        c=report['results'][name]['wall_trajectory'];ax[0,1].plot([v['t'] for v in c],[v['radius'] for v in c],label=label,ls='--' if name=='N64' else '-')
    ax[0,1].set(xlabel='Cosmic time t',ylabel='Half-height radius',title='Wall motion is measured, not prescribed');ax[0,1].legend(fontsize=9)
    im=ax[1,0].imshow((frames[-1][1][...,1]/PARAMS['fa']).T,origin='lower',extent=(-24,24,-24,24),cmap='RdBu_r',vmin=-1.1,vmax=1.1)
    ax[1,0].set(xlim=(-16,16),ylim=(-16,16),xlabel='x',ylabel='y',title='Axion angle at t = 6');fig.colorbar(im,ax=ax[1,0],label='theta = A/f')
    h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    q=np.r_[0,np.cumsum((h['Q_A'][1:]+h['Q_A'][:-1])/2*np.diff(h['t']))]
    ax[1,1].plot(h['t'],h['E_A']-h['E_A'][0],label='Axion energy change')
    ax[1,1].plot(h['t'],q,ls='--',label='Integrated transfer from S')
    ax[1,1].plot(h['t'],-(h['E_S']-h['E_S'][0]),ls=':',label='Energy lost by S')
    ax[1,1].set(xlabel='Cosmic time t',ylabel='Energy change',title='The two fields exchange energy');ax[1,1].legend(fontsize=9)
    fig.suptitle('One seeded bubble: dynamical transition field + axion\n3D CosmoLattice pilot | no expansion or friction | S_true = 1, c = 1',fontsize=13)
    fig.savefig(dest/'single_bubble.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    first=frames[0][1];images=[]
    specs=[('Transition field S','viridis',-.05,1.05),('Axion mass squared / M²','magma',0,1.15),('Axion angle theta','RdBu_r',-1.1,1.1)]
    data=lambda z:[z[...,0],z[...,0]**2,z[...,1]/PARAMS['fa']]
    for axis,values,(title,cmap,vmin,vmax) in zip(axes,data(first),specs):
        image=axis.imshow(values.T,origin='lower',extent=(-24,24,-24,24),cmap=cmap,vmin=vmin,vmax=vmax)
        images.append(image);axis.set(title=title,xlim=(-16,16),ylim=(-16,16),xlabel='x',ylabel='y');fig.colorbar(image,ax=axis,shrink=.7)
    title=fig.suptitle('Seeded bubble | t = 0 | 3D simulation, central slice')
    def update(i):
        meta,plane=frames[i]
        for image,values in zip(images,data(plane)):image.set_data(values.T)
        title.set_text(f'Seeded bubble | t = {meta["t"]:.1f} | 3D simulation, central slice')
        return images+[title]
    animation=FuncAnimation(fig,update,frames=len(frames),interval=200,blit=False)
    animation.save(dest/'single_bubble.gif',writer=PillowWriter(fps=5),dpi=95);plt.close(fig)
    from build_backend import sha
    (dest/'provenance.json').write_text(json.dumps(dict(source=str(args.source.resolve()),passed=report['passed'],scope=report['scope'],
        report_sha256=sha(args.source/'report.json'),renderer_sha256=sha(__file__)),indent=2)+'\n')
    print(dest)

if __name__=='__main__':main()

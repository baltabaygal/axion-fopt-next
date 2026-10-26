"""Visualize saved cosmological coupled fields and physically correct balances."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from model import read_state,PARAMS
from expanding_model import balances,scale,cosmic
from build_backend import sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--stochastic',action='store_true');args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'report.json').read_text());native=args.source/'main/native';h=np.genfromtxt(native/'history.csv',delimiter=',',names=True)
    meta,_=read_state(native/'initial.bin');box=meta['box'];h0=float(h['H_physical'][0]);den=box**3*PARAMS['fa']**2*PARAMS['mass']**2
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for key,label,color,style in [('K_A','Kinetic','#F0986E','-'),('G_A','Gradient','#B8A037','--'),('V_A','Potential','#5477C4','-.'),('E_A','Total','#1F2430','-')]:
        axes[0,0].plot(h['t_cosmic'],h[key]/den,color=color,ls=style,label=label)
    axes[0,0].set(xlabel='Elapsed cosmic time',ylabel=r'$\bar\rho_A/(f^2 M^2)$',title='Axion energy density');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(h['t_cosmic'],1-h['F_half'],label='All-cell S>=.5 fraction')
    axes[0,1].plot(h['t_cosmic'],h['near_true_fraction'],'--',label='Fraction |S-1|<.1')
    axes[0,1].set(xlabel='Elapsed cosmic time',ylabel='3D volume fraction',ylim=(0,1.05),title='Conversion and proximity to the minimum');axes[0,1].legend(fontsize=8)
    total,axion=balances(h);normal=max(abs(h['total'][0]),abs(h['total'][-1]),float(h['I_total'][-1]) if 'I_total' in h.dtype.names else 0,1e-30)
    axes[1,0].plot(h['t_cosmic'],total/normal,label='Total continuity residual')
    axes[1,0].plot(h['t_cosmic'],axion/normal,'--',label='Axion continuity residual')
    axes[1,0].set(xlabel='Elapsed cosmic time',ylabel='Residual / fixed energy scale',title='Expansion'+(' + insertion work' if args.stochastic else '')+' included');axes[1,0].legend(fontsize=8)
    axes[1,1].fill_between(h['t_cosmic'],h['S_min'],h['S_max'],alpha=.15,label='Global S range')
    axes[1,1].plot(h['t_cosmic'],h['S_mean'],label='Mean S')
    axes[1,1].axhline(1,color='gray',lw=.7);axes[1,1].set(xlabel='Elapsed cosmic time',ylabel='S',title='Scalar evolution in the expanding background');axes[1,1].legend(fontsize=8)
    title='Stochastic finite-radius seeding' if args.stochastic else 'One seeded bubble'
    fig.suptitle(title+f' with prescribed radiation expansion\nH0={h0:g}; a: 1 to {h["a"][-1]:.3f}; polynomial mass coupling unchanged',fontsize=13)
    fig.savefig(out/'cosmology.png',dpi=160);plt.close(fig)
    frames=[read_state(path) for path in sorted(native.glob('frame_*.bin'))]
    if len(frames)>90:frames=frames[::2]
    low=min(-.05,min(float(y[...,0].min()) for _,y in frames));high=max(1.05,max(float(y[...,0].max()) for _,y in frames))
    angle=max(1.1,max(float(abs(y[...,1]/PARAMS['fa']).max()) for _,y in frames));images=[]
    fig,axes=plt.subplots(1,3,figsize=(12,4.4),layout='constrained');data=lambda y:(y[...,0],y[...,0]**2,y[...,1]/PARAMS['fa'])
    for ax,arr,label,cmap,lo,hi in zip(axes,data(frames[0][1]),['Transition field S','Axion mass squared / M²','Axion angle theta'],['viridis','magma','RdBu_r'],[low,0,-angle],[high,max(low**2,high**2),angle]):
        im=ax.imshow(arr.T,origin='lower',extent=(-box/2,box/2,-box/2,box/2),cmap=cmap,vmin=lo,vmax=hi);images.append(im)
        ax.set(title=label,xlabel='Comoving x',ylabel='Comoving y');fig.colorbar(im,ax=ax,shrink=.7)
    heading=fig.suptitle(title+' | central slice of a 3D expanding simulation')
    def update(i):
        m,y=frames[i];u=m['t'];a=scale(u,h0);t=cosmic(u,h0)
        for im,arr in zip(images,data(y)):im.set_data(arr.T)
        heading.set_text(f'{title} | u={u:.1f}, cosmic elapsed t={t:.2f}, a={a:.3f}')
        return images+[heading]
    animation=FuncAnimation(fig,update,frames=len(frames),interval=120,blit=False);animation.save(out/'cosmology.gif',writer=PillowWriter(fps=8),dpi=90);plt.close(fig)
    (out/'provenance.json').write_text(json.dumps(dict(source=str(args.source.resolve()),report_sha256=sha(args.source/'report.json'),renderer_sha256=sha(__file__),scope='Physical energy densities; comoving plotted coordinates'),indent=2)+'\n')
    print(out)

if __name__=='__main__':main()

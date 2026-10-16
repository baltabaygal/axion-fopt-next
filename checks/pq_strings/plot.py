"""Portable scientific diagnostic figures and a phase/modulus animation."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from model import read_state, defects

BLUE='#5477C4';NEUTRAL='#464C55';ORANGE='#CC6F47'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.labelcolor':'#1F2430','text.color':'#1F2430',
                     'axes.grid':False,'savefig.facecolor':'white'})


def slice_panel(ax,state,meta,quantity):
    z=meta['n']//2
    phase=np.arctan2(state[...,1],state[...,0])
    r=np.linalg.norm(state[...,:2],axis=-1)
    data=phase if quantity=='phase' else r
    cmap='twilight_shifted' if quantity=='phase' else 'Blues'
    low,high=(-np.pi,np.pi) if quantity=='phase' else (0.,1.4)
    im=ax.imshow(data[:,:,z].T,origin='lower',extent=(0,meta['box'],0,meta['box']),
                 cmap=cmap,vmin=low,vmax=high,interpolation='nearest')
    ax.set_xlabel('Comoving x');ax.set_ylabel('Comoving y')
    diag,w=defects(state,meta['box'],meta['tau'])
    dx=meta['box']/meta['n']
    for sign,marker,color in [(1,'o','white'),(-1,'x','black')]:
        i,j=np.where(w[0][:,:,z]==sign)
        ax.scatter((i+1)*dx,(j+1)*dx,s=28,marker=marker,c=color,
                   linewidths=1.1,edgecolors='black' if marker=='o' else None)
    # Mark the theta=pi sheet only on the negative-real side; this is geometric,
    # and remains a branch-cut surface even when chi=0.
    coords=(np.arange(meta['n'])+.5)*dx
    cut=np.ma.masked_where(state[:,:,z,0].T>=0,state[:,:,z,1].T)
    if cut.count() and cut.min()<0<cut.max():
        ax.contour(coords,coords,cut,levels=[0],colors=['#F0986E'],linewidths=.8)
    ax.set_title(f"u={meta['u']:g}, tau={meta['tau']:g}")
    return im


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    report=json.loads((args.source/'summary.json').read_text())
    # Compact chart contract accompanies data and exported figures.
    contract=dict(surface='standalone scientific PNG/PDF and GIF',
                  diagnostic='periodic complex-PQ evolution in one seeded realization',
                  trends='13 physical-time snapshots; growing susceptibility versus identical massless initial data',
                  slices='fixed z=N//2 at three times, fixed cyclic phase and modulus color scales',
                  encodings='blue solid growing chi; neutral dashed massless; cyclic phase; sequential blue modulus',
                  limitations='branch-cut surface is not a physical wall when chi=0; finite box; no ensemble errors')
    (args.output/'chart_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    growing=report['runs']['network_growing'];massless=report['runs']['network_massless']
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    fig.suptitle('Small PQ network: defect evolution and energy',fontsize=16)
    fig.text(.5,.93,f"{growing['parameters']['N']}³ periodic grid; L={growing['parameters']['lSide']:g}; seed 0; n_tau={growing['parameters']['n_tau']:g}",ha='center',color='#6F768A')
    for label,run,color,style in [('Growing susceptibility',growing,BLUE,'-'),('Massless control',massless,NEUTRAL,'--')]:
        frames=run['snapshots'];u=[r['u'] for r in frames]
        axes[0,0].plot(u,[r['plaquettes'] for r in frames],color=color,ls=style,label=label,lw=1.8)
        axes[0,1].plot(u,[r['wall_links'] for r in frames],color=color,ls=style,lw=1.8)
        h=np.genfromtxt(args.source/f"network_{'growing' if run is growing else 'massless'}/native/history.csv",delimiter=',',names=True)
        axes[1,0].plot(h['u'],h['total'],color=color,ls=style,lw=1.8)
    frames=growing['snapshots']
    axes[1,1].plot([r['u'] for r in frames],[r['core_cells'] for r in frames],color=BLUE,lw=1.8,label='Core width / dx')
    axes[1,1].plot([r['u'] for r in frames],[r['wall_cells'] for r in frames],color=ORANGE,lw=1.8,label='Wall width / dx')
    axes[1,1].axhline(4,color=NEUTRAL,ls=':',lw=1.2,label='Four-cell gate')
    axes[1,1].set_yscale('log')
    titles=['String-pierced plaquettes','Links crossing theta = pi','Physical mean energy density','Characteristic widths in grid cells']
    labels=['Raw count','Raw count','rho / (f_a² omega_star²)','Width / dx']
    for ax,title,label in zip(axes.flat,titles,labels):
        ax.set_title(title,pad=10);ax.set_xlabel('Elapsed conformal time u = tau - 20');ax.set_ylabel(label)
        ax.grid(axis='y',color='#E6E8F0',lw=.6)
    axes[0,0].legend(frameon=False,fontsize=9);axes[1,1].legend(frameon=False,fontsize=8)
    fig.get_layout_engine().set(rect=(0,0,1,.90))
    for ext in ('png','pdf'):fig.savefig(args.output/f'network_diagnostics.{ext}',dpi=170)
    plt.close(fig)
    paths=sorted((args.source/'network_growing/native').glob('frame_*.bin'))
    selected=[paths[0],paths[len(paths)//2],paths[-1]]
    fig,axes=plt.subplots(2,3,figsize=(12,8),layout='constrained')
    fig.suptitle('Complex-PQ field: phase and resolved modulus',fontsize=16)
    fig.text(.5,.945,'Fixed z slice; circles/crosses mark opposite string winding; orange traces theta = pi',ha='center',fontsize=10,color='#6F768A')
    for col,path in enumerate(selected):
        meta,state=read_state(path)
        phase_image=slice_panel(axes[0,col],state,meta,'phase')
        modulus_image=slice_panel(axes[1,col],state,meta,'modulus')
    phase_bar=fig.colorbar(phase_image,ax=list(axes[0,:]),shrink=.8,label='Phase theta (rad)')
    phase_bar.set_ticks([-np.pi,0,np.pi],labels=['−π','0','π'])
    fig.colorbar(modulus_image,ax=list(axes[1,:]),shrink=.8,label='|Phi| / f_a')
    fig.get_layout_engine().set(rect=(0,0,1,.92))
    for ext in ('png','pdf'):fig.savefig(args.output/f'network_slices.{ext}',dpi=170)
    plt.close(fig)
    # GIF frames use the same z plane and color normalization throughout.
    images=[]
    for path in paths:
        meta,state=read_state(path)
        fig,axes=plt.subplots(1,2,figsize=(8,4.5),layout='constrained')
        fig.suptitle(f"PQ network, {meta['n']}³: tau={meta['tau']:g}, a={meta['a']:.2f}",fontsize=14)
        pimage=slice_panel(axes[0],state,meta,'phase');rimage=slice_panel(axes[1],state,meta,'modulus')
        fig.colorbar(pimage,ax=axes[0],shrink=.7,label='Phase (rad)');fig.colorbar(rimage,ax=axes[1],shrink=.7,label='|Phi|')
        fig.canvas.draw();images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert('RGB'))
        plt.close(fig)
    images[0].save(args.output/'network_evolution.gif',save_all=True,append_images=images[1:],duration=220,loop=0)
    # Export richer frame data behind the figures.
    (args.output/'network_frames.json').write_text(json.dumps({k:report['runs'][k]['snapshots'] for k in ('network_growing','network_massless')},indent=2)+'\n')
    print(args.output)
if __name__=='__main__':main()

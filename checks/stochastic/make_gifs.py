"""Render representative GIFs from a frozen catalogue and native field frames."""
import argparse
import hashlib
import json
import struct
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_field import ROOT,read_state,mass
from axion_fopt_next import Catalogue


def capture(fig):
    fig.canvas.draw()
    return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())


def save_gif(frames,path,duration=80):
    # Independent palettes preserve the fixed data-to-color mapping rendered
    # by Matplotlib. The final frame pauses before the loop restarts.
    frames[0].save(path,save_all=True,append_images=frames[1:],duration=[duration]*(len(frames)-1)+[600],loop=0,disposal=2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'outputs/stochastic/representative-gifs')
    parser.add_argument('--output',type=Path,default=ROOT/'checks/stochastic/results/animations')
    args=parser.parse_args();source=args.source.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((source/'animation.json').read_text()) if (source/'animation.json').exists() else dict(h=1.,beta=1.,box=3.,width=.15,origin=1.,n=32,birth_end=2.5,field_end=6.5)
    h=cfg['h'];box=cfg['box'];n=cfg['n'];origin=cfg['origin']
    z=np.load(source/'catalogue.npz');catalogue=Catalogue(box,.5,z['times'],z['positions'])
    rng=np.random.default_rng(92848);points=rng.uniform(0,box,(15000,3));arrival=catalogue.arrival_times(points)
    # A denser 2D section shows periodic sphere growth/merging clearly.
    coords=(np.arange(120)+.5)*box/120
    x,y=np.meshgrid(coords,coords,indexing='ij')
    section=catalogue.arrival_times(np.stack([x,y,np.full_like(x,box/2)],axis=-1))
    fig=plt.figure(figsize=(10,5.2),dpi=100)
    ax=fig.add_subplot(121,projection='3d');ax2=fig.add_subplot(122)
    fig.subplots_adjust(top=.78,bottom=.14,left=.04,right=.96,wspace=.2)
    times=np.linspace(0,2.7/h,55);frames=[]
    for frame,tau in enumerate(times):
        ax.clear();ax2.clear()
        converted=arrival<=tau
        # Points show the actual periodic bubble union, rather than drawing
        # overlapping transparent spheres that conceal which volume converted.
        q=points[~converted];p=points[converted]
        if len(q):ax.scatter(*q.T,s=1,color='#b8c1ca',alpha=.035,depthshade=False,rasterized=True)
        if len(p):ax.scatter(*p.T,s=3,color='#19a7a0',alpha=.5,depthshade=True,rasterized=True)
        for a in (0,box):
            for b in (0,box):
                ax.plot([0,box],[a,a],[b,b],color='#485461',lw=.8)
                ax.plot([a,a],[0,box],[b,b],color='#485461',lw=.8)
                ax.plot([a,a],[b,b],[0,box],color='#485461',lw=.8)
        ax.set(xlim=(0,box),ylim=(0,box),zlim=(0,box),xlabel='x',ylabel='y',zlabel='z')
        ax.set_box_aspect((1,1,1));ax.view_init(elev=22,azim=-55);ax.set_title(f'Periodic comoving cube, L={box:g}')
        from matplotlib.colors import ListedColormap
        ax2.imshow((section<=tau).T,origin='lower',extent=(0,box,0,box),vmin=0,vmax=1,
                   cmap=ListedColormap(['#e8edf2','#19a7a0']),interpolation='nearest')
        ax2.set(xlabel='x',ylabel='y',title=f'Slice at z={box/2:g}');ax2.set_aspect('equal')
        fig.suptitle('Stochastic bubble nucleation and growth',fontsize=17,y=.98)
        count=int(np.sum(catalogue.times<=tau))
        fig.text(.5,.885,f'Conformal time tau={tau:.2f}   |   Cosmic time t={h*tau**2/2:.2f}   |   Bubbles born: {count}',ha='center',fontsize=11)
        fig.text(.5,.835,f'True vacuum: teal   ·   False vacuum: light grey   ·   False fraction ≈ {np.mean(~converted):.3f}',ha='center',fontsize=10)
        fig.text(.5,.03,f'H*={h:g} · beta={cfg["beta"]:g} · v=0.5 · 100,000 tracers · seed 0 · comoving axes',ha='center',fontsize=10)
        frames.append(capture(fig))
        # Remove figure text before rendering the next frame.
        for text in list(fig.texts):text.remove()
        if frame%15==0:print(f'Bubble frame {frame+1}/{len(times)}',flush=True)
    save_gif(frames,out/'bubble_nucleation.gif',100)
    frames[28].save(out/'bubble_preview.png');plt.close(fig);del frames
    files=sorted((source/'evolution').glob('frame_*.bin'))
    history=np.genfromtxt(source/'evolution/history.csv',delimiter=',',names=True)
    fig,(ax,trace)=plt.subplots(1,2,figsize=(10,4.9),dpi=100,gridspec_kw={'width_ratios':[1,1.25]})
    fig.subplots_adjust(top=.8,bottom=.17,left=.07,right=.95,wspace=.6)
    image=ax.imshow(np.zeros((n,n)),origin='lower',extent=(0,box,0,box),vmin=-1,vmax=1,cmap='RdBu_r',interpolation='bilinear')
    fig.colorbar(image,ax=ax,label='theta',fraction=.046,pad=.04)
    ax.set(xlabel='x (comoving)',ylabel='y (comoving)',title=f'Field slice at z={(n//2+.5)*box/n:.3f}')
    t=h*history['tau']**2/2
    trace.plot(t,history['theta_mean'],color='#263f60',lw=1.7,label='Spatial mean theta')
    # A spatial average can cancel while the local field keeps oscillating.
    # Read one fixed site's actual native values without loading whole frames.
    local_times=[];local_values=[]
    for file in files:
        with file.open('rb') as stream:
            header=stream.read(36)
            dimension,_,u,_=struct.unpack('<Iddd',header[8:])
            if dimension!=n:raise AssertionError('Frame dimension mismatch')
            index=(n//4*n+n//4)*n+n//2
            stream.seek(36+32*index)
            local_values.append(float(np.fromfile(stream,dtype='<f8',count=1)[0]))
            local_times.append(h*(origin+u)**2/2)
    trace.plot(local_times,local_values,color='#9467bd',lw=1.,alpha=.85,label='Local theta (fixed point)')
    trace.set(xlim=(float(t[0]),float(t[-1])),ylim=(-1.05,1.05),xlabel='Cosmic time t=H* tau²/2',ylabel='Mean theta')
    trace.grid(alpha=.18);trace.axhline(0,color='#9aa5b1',lw=.8)
    line=trace.axvline(t[0],color='#e49b32',lw=1.5);dot,=trace.plot([t[0]],[1.],marker='o',color='#e49b32')
    trace.legend(loc='upper right',fontsize=9)
    fig.suptitle('Axion field oscillations after stochastic nucleation',fontsize=16,y=.97)
    label=fig.text(.5,.86,'',ha='center',fontsize=11)
    fig.text(.5,.025,f'H*={h:g} · mt=1 · cosine · N={n} · theta(tau={origin:g})=1 · births end at tau={cfg["birth_end"]:g}',ha='center',fontsize=9)
    frames=[];means=[]
    for i,file in enumerate(files):
        meta,state=read_state(file);tau=origin+meta['time'];cosmic=h*tau**2/2
        if abs(meta['scale']-h*tau)>1e-11:raise AssertionError('Scale factor mismatch')
        expected=.5*(1+np.tanh((tau-state[...,2])*.5/cfg['width']))
        np.testing.assert_allclose(state[...,3],expected,atol=3e-14,rtol=0)
        image.set_data(state[:,:,n//2,0].T);line.set_xdata([cosmic,cosmic]);dot.set_data([cosmic],[state[...,0].mean()])
        label.set_text(f'Conformal time tau={tau:.2f}   |   Cosmic time t={cosmic:.2f}   |   Scale factor a={h*tau:.2f}')
        frames.append(capture(fig));means.append(float(state[...,0].mean()))
        if i%30==0:print(f'Field frame {i+1}/{len(files)}',flush=True)
    # Require genuinely oscillating simulated data, not just decreasing theta.
    crossings=int(np.count_nonzero(np.diff(np.sign(means))))
    if crossings<3:raise AssertionError('Too few resolved oscillations')
    save_gif(frames,out/'axion_oscillations.gif',75)
    frames[40].save(out/'field_preview.png');plt.close(fig)
    # The extension must reproduce the earlier validated interval exactly.
    same_early=None
    if h==1. and n==32:
        _,check=read_state(source/'evolution/frame_0000300.bin')
        _,old=read_state(ROOT/'outputs/stochastic/field-reference/N32_dt0.005/final.bin')
        np.testing.assert_array_equal(check,old);same_early=True
    metadata=dict(bubble_frames=len(times),field_frames=len(files),mean_zero_crossings=crossings,
                  catalogue_bubbles=len(catalogue.times),source_directory=str(source),
                  physics=cfg,
                  limits='Representative visualization. Frozen birth catalogue, continued bubble growth, prescribed geometry, no backreaction. Cube uses independent probes; spatial continuum accuracy is not certified.',
                  same_early_field_exact=same_early,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  executable_sha256=hashlib.sha256((ROOT/'build/cosmolattice/axion_fopt_evolve').read_bytes()).hexdigest(),
                  input_sha256=hashlib.sha256((source/'input.json').read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':main()

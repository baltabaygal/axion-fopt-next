"""Show below-threshold cells throughout the volume, not just the central slice."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import read_state,PARAMS
from expanding_model import physical_terms,scale,cosmic

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    m,y=read_state(args.source/'main/native/final.bin');h0=json.loads((args.source/'main/cosmology.json').read_text())['H0']
    n=m['n'];a=float(scale(m['t'],h0));box=m['box'];dx=box/n
    s=y[...,0];mask=s<.5;minimum=np.unravel_index(np.argmin(s),s.shape);k=minimum[2]
    depth=mask.sum(axis=2)*a*dx;values=physical_terms(y,dx,a);ea=values['K_A']+values['G_A']+values['V_A']
    extent=(-a*box/2,a*box/2,-a*box/2,a*box/2)
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    data=[s[:,:,k],depth,ea[:,:,k]/(PARAMS['fa']**2*PARAMS['mass']**2)]
    titles=[f'S through the global minimum\nphysical z={(k+.5)*a*dx-a*box/2:.2f}',
        'Below-threshold depth through the full box\nIntegral of [S<0.5] along physical z',r'Axion density on the minimum plane'+'\n'+r'$\rho_A/(f^2M^2)$']
    for ax,v,title,cmap in zip(axes,data,titles,['viridis','magma','cividis']):
        im=ax.imshow(v.T,origin='lower',extent=extent,cmap=cmap);ax.set(title=title,xlabel='Physical x',ylabel='Physical y');fig.colorbar(im,ax=ax,shrink=.65)
    if s[:,:,k].min()<.5<s[:,:,k].max():
        coord=(np.arange(n)+.5)*a*dx-a*box/2
        axes[0].contour(coord,coord,s[:,:,k].T,levels=[.5],colors='white',linewidths=.6)
    fig.suptitle(f'Final 3D remnant diagnostic | u={m["t"]:.2f}, cosmic elapsed t={float(cosmic(m["t"],h0)):.2f}, a={a:.3f}\n'
        f'{int(mask.sum())} of {n**3} cells below S=0.5; white contour marks the threshold; periodic box',fontsize=12)
    fig.savefig(args.output,dpi=170);plt.close(fig)

if __name__=='__main__':main()

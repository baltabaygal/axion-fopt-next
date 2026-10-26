"""Continuum-defined symmetric standing fields in a biased planar slab.

Smoothly tapered sine modes are initial data, not cavity eigenmodes. The
pressure normalization is a linear, untapered phase-average estimate only.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from paper_reference import Parameters,terms,head_on_reflection
from planar_stencils import gradient_density

@dataclass(frozen=True)
class SlabSeed:
    box:float=512.
    length:float=256.
    width:float=1.
    speed:float=.05
    taper:float=16.
    pressure_ratio:float=1.


def prepare(n,seed=SlabSeed(),p=Parameters(dv=.0002)):
    if (not isinstance(n,int) or not 128<=n<=16384 or
        not np.isfinite(list(asdict(seed).values())).all() or
        not 0<seed.length<seed.box or seed.width<=0 or not 0<=seed.speed<1 or
        not 0<seed.taper<seed.length/4 or seed.pressure_ratio<0 or p.dv<=0 or
        abs(p.lam-8/seed.width**2)>1e-12):raise ValueError('Invalid distributed slab seed')
    indices=np.array([32,36,40]);phases=np.array([0.,.7,1.4])
    k=np.pi*indices/seed.length;omega=np.sqrt(k*k+p.mass_false**2)
    if k.max()>=np.pi*n/seed.box:raise ValueError('Unresolved seed modes')
    amplitude=np.sqrt(4*p.dv*seed.pressure_ratio/np.sum(k*k))
    x=(np.arange(n)+.5)*seed.box/n;left=(seed.box-seed.length)/2;right=left+seed.length
    gamma=1/np.sqrt(1-seed.speed**2);s=np.zeros((n,4));s[:,0]=1.
    for image in (-seed.box,0.,seed.box):
        for location,sign,velocity in ((left,-1.,seed.speed),(right,1.,-seed.speed)):
            z=gamma*(x-location-image)/seed.width;t=np.tanh(z)
            s[:,0]+=.5*sign*t;s[:,2]+=-velocity*.5*sign*gamma/seed.width*(1-t*t)
    distance=np.minimum(x-left,right-x)
    q=np.clip((distance-8*seed.width)/seed.taper,0,1)
    window=q**3*(10-15*q+6*q*q)
    for kk,ww,phase in zip(k,omega,phases):
        shape=window*np.sin(kk*(x-left))
        s[:,1]+=amplitude*np.cos(phase)*shape
        s[:,3]+=-amplitude*ww*np.sin(phase)*shape
    return x,s,dict(indices=indices.tolist(),phases=phases.tolist(),k=k.tolist(),omega=omega.tolist(),
        amplitude_per_mode=float(amplitude),linear_phase_average_pressure=float(amplitude**2*np.sum(k*k)/4),
        taper_policy='Quintic rise from zero at distance 8 widths to unity 16 units farther inside each wall',
        interpretation='Coherent standing-field superposition; not thermal, isotropic or exact eigenmodes')


def audit(n,seed=SlabSeed(),p=Parameters(dv=.0002)):
    x,s,m=prepare(n,seed,p);dx=seed.box/n;v=terms(s[:,0],s[:,1],p)[1]
    energy=.5*s[:,3]**2+v+gradient_density(s[:,1],dx,4)
    wallfreq=[head_on_reflection(w,p.mass_false,seed.speed)['wall_frame_frequency'] for w in m['omega']]
    kfft=2*np.pi*np.fft.fftfreq(n,d=dx);wfft=np.sqrt(kfft*kfft+p.mass_false**2)
    spectral=.5*(abs(np.fft.fft(s[:,3]))**2+wfft*wfft*abs(np.fft.fft(s[:,1]))**2)
    shifted=(wfft+seed.speed*abs(kfft))/np.sqrt(1-seed.speed**2)
    mt=p.mass_false*np.exp(p.kappa/2)
    above=float(spectral[shifted>=mt].sum()/spectral.sum()) if spectral.sum() else 0.
    return dict(n=n,dx=dx,seed=asdict(seed),parameters=asdict(p),modes=m,
        initial_fv_fourier_energy_fraction_above_boosted_threshold=above,
        axion_energy_per_area=float(np.sum(energy)*dx),
        tv_energy_per_area=float(np.sum(energy[s[:,0]>.9])*dx),
        axion_angle_max=float(np.max(abs(s[:,1]))/p.fa),
        reflection_symmetry_error=float(np.max(abs(s[:,0]-s[::-1,0]))),
        axion_antisymmetry_error=float(np.max(abs(s[:,[1,3]]+s[::-1,[1,3]]))),
        wall_frame_carrier_frequencies=wallfreq,
        carrier_penetration_depths=[float(1/np.sqrt((p.mass_false*np.exp(p.kappa/2))**2-w*w)) for w in wallfreq],
        caveat='Taper broadens spectrum; carrier thresholds do not guarantee confinement of the full field')


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[audit(n) for n in (2048,4096,8192)]
    result=dict(scope='Initial data only, no solver evolution',rows=rows,solver_runs=0,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()

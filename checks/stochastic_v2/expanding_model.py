"""Independent conformal-time reference and inputs for prescribed radiation era."""
import os
from pathlib import Path
import struct
import subprocess
import time
import numpy as np
from model import PARAMS,terms,rhs,read_state

def scale(u,h0):return 1+h0*np.asarray(u)
def cosmic(u,h0):return np.asarray(u)+h0*np.asarray(u)**2/2

def expanding_rhs(u,flat,n,box,h0):
    y=flat.reshape(n,n,n,4);a=scale(u,h0);dx=box/n
    ordinary=rhs(u,flat,n,box).reshape(y.shape)
    fields=y[...,:2]
    lap=sum(np.roll(fields,1,i)+np.roll(fields,-1,i)-2*fields for i in range(3))/dx**2
    force=lap-ordinary[...,2:]
    return np.concatenate((y[...,2:]/a**2,a**2*lap-a**4*force),axis=-1).ravel()

def physical_terms(y,dx,a):
    values=terms(y,dx)
    return {key:v/(a**6 if key.startswith('K_') else a**2 if key.startswith('G_') else 1) for key,v in values.items()}

def cumulative(y,x):return np.r_[0,np.cumsum((y[1:]+y[:-1])/2*np.diff(x))]

def balances(h):
    u=h['u'];ka=h['K_A'];ks=h['K_S'];ga=h['G_A'];gs=h['G_S'];hc=h['H_conformal']
    total_integral=cumulative(-hc*(6*(ka+ks)+2*(ga+gs)),u)
    axion_integral=cumulative(h['Q_A']-hc*(6*ka+2*ga),u)
    jumps_total=h['I_total']-h['I_total'][0] if 'I_total' in h.dtype.names else 0.
    jumps_axion=h['I_A']-h['I_A'][0] if 'I_A' in h.dtype.names else 0.
    return h['total']-h['total'][0]-total_integral-jumps_total,h['E_A']-h['E_A'][0]-axion_integral-jumps_axion

def invoke(executable,directory,y,box,dt,end,h0,cadence=.04,frames=.2,extra=None,timeout=180):
    dest=Path(directory).resolve();dest.mkdir(parents=True,exist_ok=False);n=len(y)
    with (dest/'input.bin').open('wb') as f:
        f.write(struct.pack('<8sId',b'PKIN0001',n,box));f.write(np.asarray(y,dtype='<f8').tobytes())
    cfg=dict(N=n,lSide=box,dt=dt,tMax=end,t0=0,a0=1,expansion=str(h0>0).lower(),fixedBackground='true',H0=h0,omegaEoS=1/3,
        withGWs='false',evolver='PV2',delta_v=PARAMS['dv'],mass_tv=PARAMS['mass'],axion_fa=PARAMS['fa'],
        output_every=max(1,round(cadence/dt)),snapshot_every=max(1,round(frames/dt)),
        initial_state=str(dest/'input.bin'),result_dir=str(dest/'native'),PS_type=1)
    cfg['lambda']=PARAMS['lam'];cfg.update(extra or {})
    (dest/'input.in').write_text('\n'.join(f'{k} = {v}' for k,v in cfg.items())+'\n')
    (dest/'cosmology.json').write_text(__import__('json').dumps(dict(H0=h0,time='elapsed conformal u; t_cosmic=u+H0*u^2/2',momentum='P=a^2 X_prime=a^3 X_dot',energy='physical densities integrated over comoving coordinate volume'),indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    started=time.monotonic();p=subprocess.run([str(Path(executable).resolve()),'input='+str(dest/'input.in')],env=env,capture_output=True,text=True,timeout=timeout)
    (dest/'run.log').write_text(p.stdout+p.stderr)
    if p.returncode:raise RuntimeError(p.stderr[-3000:])
    return time.monotonic()-started

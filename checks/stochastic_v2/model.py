"""Independent expanded-potential reference, initial data and binary interchange."""
import os
from pathlib import Path
import struct
import subprocess
import time
import numpy as np

PARAMS=dict(lam=2.,dv=.06,mass=.5,fa=.2)

def terms(y,dx,p=PARAMS):
    s,a,ps,pa=np.moveaxis(y,-1,0)
    u=p['lam']/4*(s*s-2*s**3+s**4)+p['dv']*(1-4*s**3+3*s**4)
    va=p['fa']**2*p['mass']**2*s*s*(1-np.cos(a/p['fa']))
    gs=sum((np.roll(s,-1,i)-s)**2 for i in range(3))/(2*dx*dx)
    ga=sum((np.roll(a,-1,i)-a)**2 for i in range(3))/(2*dx*dx)
    return dict(K_S=ps*ps/2,G_S=gs,U=u,K_A=pa*pa/2,G_A=ga,V_A=va)

def rhs(t,flat,n,box,p=PARAMS):
    y=flat.reshape(n,n,n,4);s,a,ps,pa=np.moveaxis(y,-1,0);dx=box/n
    fs=p['lam']/4*(2*s-6*s*s+4*s**3)+p['dv']*(-12*s*s+12*s**3)+2*p['fa']**2*p['mass']**2*s*(1-np.cos(a/p['fa']))
    fa=p['fa']*p['mass']**2*s*s*np.sin(a/p['fa'])
    lap=lambda x:sum(np.roll(x,1,i)+np.roll(x,-1,i)-2*x for i in range(3))/dx**2
    return np.stack((ps,pa,lap(s)-fs,lap(a)-fa),-1).ravel()

def read_state(path):
    with Path(path).open('rb') as stream:
        magic,n,box,t=struct.unpack('<8sIdd',stream.read(28))
        if magic not in (b'PKST0001',b'PKSL0001'):raise ValueError('Bad state header')
        shape=(n,n,n,4) if magic==b'PKST0001' else (n,n,4)
        y=np.frombuffer(stream.read(),dtype='<f8').copy().reshape(shape)
    return dict(n=n,box=box,t=t),y

def seed(n,box,theta=1.):
    x=(np.arange(n)+.5)*box/n-box/2
    r=np.sqrt(x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2)
    y=np.zeros((n,n,n,4));y[...,0]=.5*(1-np.tanh((r-8)/2));y[...,1]=PARAMS['fa']*theta
    return y

def run(executable,directory,y,box,dt,end,cadence=.08,frames=.4,params=PARAMS,timeout=55):
    directory=Path(directory).resolve();directory.mkdir(parents=True,exist_ok=False);n=len(y)
    with (directory/'input.bin').open('wb') as stream:
        stream.write(struct.pack('<8sId',b'PKIN0001',n,box));stream.write(np.asarray(y,dtype='<f8').tobytes())
    cfg=dict(N=n,lSide=box,dt=dt,tMax=end,t0=0,a0=1,expansion='false',fixedBackground='true',
        withGWs='false',evolver='PV2',lambda_=params['lam'],delta_v=params['dv'],mass_tv=params['mass'],axion_fa=params['fa'],
        output_every=max(1,round(cadence/dt)),snapshot_every=max(1,round(frames/dt)),
        initial_state=str(directory/'input.bin'),result_dir=str(directory/'native'),kIR=2*np.pi/box,
        omegaEoS=0,PS_type=1)
    cfg['lambda']=cfg.pop('lambda_')
    (directory/'input.in').write_text('\n'.join(f'{k} = {v}' for k,v in cfg.items())+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='false',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    start=time.monotonic()
    result=subprocess.run([str(Path(executable).resolve()),'input='+str(directory/'input.in')],env=env,capture_output=True,text=True,timeout=timeout)
    (directory/'run.log').write_text(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError(result.stderr[-3000:])
    return time.monotonic()-start

"""Independent planar wall geometry and exact semi-discrete energy accounting."""
import numpy as np
from paper_reference import terms


def wall_geometry(state,box):
    state=np.asarray(state);n=len(state);dx=box/n
    if state.shape!=(n,4) or not np.isfinite(state).all() or not np.isfinite(box) or box<=0:
        raise ValueError('Require finite line state and positive box')
    phi=state[:,0]
    def crossing(level):
        following=np.roll(phi,-1)
        cells=np.flatnonzero((phi<=level)&(following>level))
        if len(cells)!=1:raise ValueError('Ambiguous rising wall level crossing')
        i=int(cells[0]);j=(i+1)%n;fraction=(level-phi[i])/(phi[j]-phi[i])
        position=((i+.5+fraction)*dx)%box
        momentum=(1-fraction)*state[i,2]+fraction*state[j,2]
        return position,-momentum/((phi[j]-phi[i])/dx)
    centre,velocity=crossing(.5);low,_=crossing(.1);high,_=crossing(.9)
    nearest=lambda d:(d+box/2)%box-box/2
    lower=nearest(low-centre);upper=nearest(high-centre)
    if not lower<0<upper:raise ValueError('Invalid wall width ordering')
    return dict(position=centre,phase_velocity=velocity,width_10_90=upper-lower,
                width_cells=(upper-lower)/dx)


def window_budget(state,box,weight,weight_dot,p):
    """Signed flux includes boundary advection; Q positive into axion sector."""
    state=np.asarray(state);n=len(state);dx=box/n
    weight=np.asarray(weight);weight_dot=np.asarray(weight_dot)
    if state.shape!=(n,4) or weight.shape!=(n,) or weight_dot.shape!=(n,) or not np.isfinite(state).all() or not np.isfinite(weight).all() or not np.isfinite(weight_dot).all():
        raise ValueError('Invalid line window data')
    phi,a,pp,pa=state.T;u,va,_,_=terms(phi,a,p)
    values=[]
    for field,momentum,potential in ((phi,pp,u),(a,pa,va)):
        difference=np.roll(field,-1)-field
        density=.5*momentum**2+potential+.25*(difference**2+np.roll(difference,1)**2)/dx**2
        current=-difference*(momentum+np.roll(momentum,-1))/(2*dx)
        lab_flux=float(np.sum(current*(weight-np.roll(weight,-1))))
        advection=float(np.sum(weight_dot*density)*dx)
        values.append((float(np.sum(weight*density)*dx),lab_flux-advection,lab_flux,advection))
    return dict(energy=values[0][0]+values[1][0],phi_energy=values[0][0],axion_energy=values[1][0],
        phi_flux=values[0][1],axion_flux=values[1][1],
        phi_lab_flux=values[0][2],axion_lab_flux=values[1][2],
        phi_advection=values[0][3],axion_advection=values[1][3],
        exchange=float(np.sum(weight*p.kappa*va*pp)*dx))

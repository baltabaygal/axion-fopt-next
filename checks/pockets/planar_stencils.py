"""Independent periodic spatial stencils and matching Hamiltonian accounting.

Planning/reference module only; the compiled line driver still uses order 2.
Fourth-order signed bonds give nonnegative *total* gradient energy, not a
pointwise positive split density. Do not interpret negative bond terms as
physical negative energy or silently reuse nearest-neighbour flux formulas.
"""
import numpy as np
from paper_reference import terms


def bonds(order):
    if order==2:return ((1,1.),)
    if order==4:return ((1,4./3.),(2,-1./12.))
    raise ValueError('Spatial order must be 2 or 4')


def laplacian(field,dx,order=2):
    field=np.asarray(field)
    if field.ndim!=1 or len(field)<5 or not np.isfinite(field).all() or not np.isfinite(dx) or dx<=0:
        raise ValueError('Require a finite periodic line and positive dx')
    return sum(c*(np.roll(field,-hop)+np.roll(field,hop)-2*field)/dx**2 for hop,c in bonds(order))


def gradient_density(field,dx,order=2):
    """Endpoint-shared bond Hamiltonian, including signed order-4 hop-2 term."""
    return sum(c*((np.roll(field,-hop)-field)**2+(field-np.roll(field,hop))**2)/(4*dx**2) for hop,c in bonds(order))


def symbol(k,dx,order=2):
    """Eigenvalue of minus Laplacian; nonnegative for both supported orders."""
    if not np.isfinite(dx) or dx<=0:raise ValueError('Positive finite dx required')
    return sum(4*c*np.sin(hop*np.asarray(k)*dx/2)**2/dx**2 for hop,c in bonds(order))


def window_budget(state,box,weight,weight_dot,p,order=2):
    state=np.asarray(state);weight=np.asarray(weight);weight_dot=np.asarray(weight_dot)
    n=len(state);dx=box/n
    if state.shape!=(n,4) or weight.shape!=(n,) or weight_dot.shape!=(n,) or not np.isfinite(box) or box<=0 or not np.isfinite(state).all() or not np.isfinite(weight).all() or not np.isfinite(weight_dot).all():
        raise ValueError('Invalid spatial-reference window data')
    phi,a,pp,pa=state.T;u,va,_,_=terms(phi,a,p);result=[]
    for field,momentum,potential in ((phi,pp,u),(a,pa,va)):
        density=.5*momentum**2+potential+gradient_density(field,dx,order)
        flux=0.
        for hop,c in bonds(order):
            current=-c*(np.roll(field,-hop)-field)*(momentum+np.roll(momentum,-hop))/(2*dx)
            flux+=float(np.sum(current*(weight-np.roll(weight,-hop))))
        advection=float(np.sum(weight_dot*density)*dx)
        result.append((float(np.sum(weight*density)*dx),flux-advection,flux,advection))
    return dict(energy=result[0][0]+result[1][0],phi_energy=result[0][0],axion_energy=result[1][0],
        phi_flux=result[0][1],axion_flux=result[1][1],phi_lab_flux=result[0][2],axion_lab_flux=result[1][2],
        phi_advection=result[0][3],axion_advection=result[1][3],exchange=float(np.sum(weight*p.kappa*va*pp)*dx))

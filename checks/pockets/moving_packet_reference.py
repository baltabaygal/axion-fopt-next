"""Lorentz-boosted weak-axion field reference on a degenerate moving kink."""
import math
import numpy as np
from planar_initial import PlanarSeed,modes
from paper_reference import Parameters,head_on_reflection
from formation_scattering_reference import scatter


def reference_modes(seed=PlanarSeed(speed=.3),p=Parameters(dv=0.),quadrature_factor=1):
    indices,k,omega,weights=modes(seed,p)
    if not isinstance(quadrature_factor,int) or quadrature_factor<1:raise ValueError('Positive integer quadrature factor required')
    if quadrature_factor>1:
        k=np.linspace(k[0],k[-1],(len(k)-1)*quadrature_factor+1)
        k0=math.sqrt(seed.omega**2-p.mass_false**2)
        weights=np.exp(-.5*(seed.envelope_width*(k-k0))**2)
        weights[[0,-1]]*=.5;weights/=weights.sum()
        omega=np.sqrt(k*k+p.mass_false**2)
    gamma=1/math.sqrt(1-seed.speed**2)
    rows=[]
    for kk,ww,weight in zip(k,omega,weights):
        rest_omega=gamma*(ww+seed.speed*kk)
        row=scatter(float(rest_omega),mass_false=p.mass_false,kappa=p.kappa,width=seed.width,samples=513)
        row.update(k=float(kk),lab_omega=float(ww),weight=float(weight),
            rest_k=float(gamma*(kk+seed.speed*ww)),quadrature_period=seed.box*quadrature_factor,
            reflected_lab_omega=head_on_reflection(float(ww),p.mass_false,seed.speed)['reflected_lab_frequency'])
        rows.append(row)
    return rows


def packet_fields(x,t,rows,seed=PlanarSeed(speed=.3),wall_initial=256.):
    gamma=1/math.sqrt(1-seed.speed**2);z=np.asarray(x)-wall_initial
    xi=gamma*(z+seed.speed*t);tau=gamma*(t+seed.speed*z)
    a=np.zeros_like(z,dtype=float);pa=np.zeros_like(a);gradient=np.zeros_like(a)
    for row in rows:
        kk=row['k'];om=row['lab_omega'];rest=row['omega'];rk=row['rest_k']
        reflection=complex(row['reflection_real'],row['reflection_imag'])
        kt=np.sqrt(complex(rest**2-row['mass_true']**2))
        psi=np.empty_like(xi,dtype=complex);derivative=np.empty_like(psi)
        left=xi<row['x'][0];right=xi>row['x'][-1];middle=~(left|right)
        psi[left]=np.exp(1j*rk*xi[left])+reflection*np.exp(-1j*rk*xi[left])
        derivative[left]=1j*rk*(np.exp(1j*rk*xi[left])-reflection*np.exp(-1j*rk*xi[left]))
        psi[right]=row['wave'][-1]*np.exp(1j*kt*(xi[right]-row['x'][-1]));derivative[right]=1j*kt*psi[right]
        psi[middle]=np.interp(xi[middle],row['x'],row['wave']);derivative[middle]=np.interp(xi[middle],row['x'],row['derivative'])
        c=seed.amplitude*row['weight']*np.exp(1j*kk*(wall_initial-seed.centre)-1j*rest*tau)
        a+=(c*psi).real
        pa+=(c*gamma*(seed.speed*derivative-1j*rest*psi)).real
        gradient+=(c*gamma*(derivative-1j*seed.speed*rest*psi)).real
    return a,pa,gradient,.5*(1+np.tanh(xi/seed.width))


def left_asymptotic_fields(x,t,predictions,seed=PlanarSeed(speed=.3),wall_initial=256.):
    """Exact asymptotic modes; caller must remain well left of the moving wall."""
    x=np.asarray(x);a=np.zeros_like(x,dtype=float);pa=np.zeros_like(a)
    for r in predictions:
        k=r['k'];om=r['lab_omega'];out=r['reflected_lab_omega']
        kout=np.sqrt(out*out-.1**2)
        reflection=complex(r['reflection_real'],r['reflection_imag'])
        incoming=seed.amplitude*r['weight']*np.exp(1j*k*(x-seed.centre)-1j*om*t)
        returned=seed.amplitude*r['weight']*reflection*np.exp(1j*k*(wall_initial-seed.centre)-1j*kout*(x-wall_initial)-1j*out*t)
        a+=(incoming+returned).real;pa+=(-1j*om*incoming-1j*out*returned).real
    return a,pa

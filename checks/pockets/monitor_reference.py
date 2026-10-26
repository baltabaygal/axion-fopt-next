"""Independent semi-discrete energy flux through a fixed spherical monitor.

Bond gradient energy is split equally between its endpoints. The face
current -grad(field)*(momentum_left+momentum_right)/2 then obeys an exact
semi-discrete continuity identity. This is net flux, including any incoming
radiation, not irreversible escaped energy. Static canonical momenta only.
"""
import numpy as np

from model import dilaton_potential, radius


def monitor(state, box, monitor_radius, lam, dv, mass2, mass_power=2, vacuum_stiffening=0.,formation=None):
    n=len(state);dx=box/n
    if state.shape!=(n,n,n,4) or not 0<monitor_radius<box/2 or mass_power not in (2,4):
        raise ValueError('Require cubic four-component state and a monitor inside half the box')
    mask=radius(n,box)<monitor_radius
    phi,a,pphi,pa=(state[...,j] for j in range(4))
    phi_energy=.5*pphi**2+dilaton_potential(phi,lam,dv,vacuum_stiffening)
    a_energy=.5*pa**2+.5*mass2*phi**mass_power*a*a
    if formation is not None:
        from paper_reference import terms
        if formation.lam!=lam or formation.dv!=dv or vacuum_stiffening!=0.:
            raise ValueError('Inconsistent formation vacuum parameters')
        u,va,_,_=terms(phi,a,formation)
        phi_energy=.5*pphi**2+u;a_energy=.5*pa**2+va
    fluxes=[]
    for field,momentum,density in ((phi,pphi,phi_energy),(a,pa,a_energy)):
        flux=0.
        for axis in range(3):
            following=np.roll(field,-1,axis);previous=np.roll(field,1,axis)
            density+=((following-field)**2+(field-previous)**2)/(4*dx**2)
            direction=mask.astype(np.int8)-np.roll(mask,-1,axis).astype(np.int8)
            current=-(following-field)*(momentum+np.roll(momentum,-1,axis))/(2*dx)
            flux+=float(np.sum(current*direction)*dx**2)
        fluxes.append(flux)
    energy_phi=float(np.sum(phi_energy[mask])*dx**3)
    energy_a=float(np.sum(a_energy[mask])*dx**3)
    exchange=float(np.sum((.5*mass_power*mass2*phi**(mass_power-1)*a*a*pphi)[mask])*dx**3)
    if formation is not None:
        exchange=float(np.sum((formation.kappa*va*pphi)[mask])*dx**3)
    return dict(monitor_energy=energy_phi+energy_a,monitor_phi_energy=energy_phi,
                monitor_axion_energy=energy_a,monitor_phi_flux=fluxes[0],
                monitor_axion_flux=fluxes[1],monitor_exchange=exchange)


def continuity_summary(h, prefix='monitor'):
    """Trapezoidal output-cadence audit, with a two-times-coarser comparison."""
    if prefix not in ('monitor','tracking'):
        raise ValueError('Unknown energy-current monitor')
    if prefix=='tracking':
        names=('energy','phi_energy','axion_energy','phi_flux','axion_flux','exchange')
        if 'tracking_energy' not in h.dtype.names:
            return None
        selected=np.empty(h.shape,dtype=[('u',float)]+[('monitor_'+name,float) for name in names])
        selected['u']=h['u']
        for name in names:
            selected['monitor_'+name]=h['tracking_'+name]
        result=continuity_summary(selected)
        result['interpretation']='Moving-window net current includes boundary advection and incoming waves; not irreversible escape. Cadence comparison is not time or spatial refinement.'
        return result
    if 'monitor_energy' not in h.dtype.names:
        return None
    t=h['u'];scale=abs(h['monitor_energy'][0])
    def integrate(values):
        return np.concatenate(([0.],np.cumsum(.5*(values[1:]+values[:-1])*np.diff(t))))
    phi_flux=integrate(h['monitor_phi_flux']);a_flux=integrate(h['monitor_axion_flux'])
    exchange=integrate(h['monitor_exchange'])
    residuals=dict(total=(h['monitor_energy']-h['monitor_energy'][0]+phi_flux+a_flux)/scale,
                   phi=(h['monitor_phi_energy']-h['monitor_phi_energy'][0]+phi_flux+exchange)/scale,
                   axion=(h['monitor_axion_energy']-h['monitor_axion_energy'][0]+a_flux-exchange)/scale)
    indices=np.arange(0,len(h),2)
    if indices[-1]!=len(h)-1:
        indices=np.append(indices,len(h)-1)
    cadence={key:float(abs(np.trapezoid(h[key],t)-np.trapezoid(h[key][indices],t[indices]))/scale)
             for key in ('monitor_phi_flux','monitor_axion_flux','monitor_exchange')}
    return dict(energy_scale=float(scale),net_outward_phi_energy=float(phi_flux[-1]),
                net_outward_axion_energy=float(a_flux[-1]),net_energy_into_axions=float(exchange[-1]),
                inward_phi_energy=float(np.trapezoid(np.maximum(-h['monitor_phi_flux'],0),t)),
                inward_axion_energy=float(np.trapezoid(np.maximum(-h['monitor_axion_flux'],0),t)),
                max_continuity_residual={key:float(np.max(np.abs(value))) for key,value in residuals.items()},
                final_continuity_residual={key:float(value[-1]) for key,value in residuals.items()},
                coarser_cadence_difference=cadence,
                interpretation='Net flux includes incoming waves. Cadence comparison is not time or spatial refinement.')

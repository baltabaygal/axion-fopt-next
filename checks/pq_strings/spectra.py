"""Finite-volume axion-current spectra inspired by arXiv:2609.20801 Eqs.7,21,24.

Currents use the modulus-weighted Cartesian projection. Spatial currents live
on forward-link midpoints; their FFTs are aligned by exp(-ik_i dx/2) before
projection onto K_i=2*sin(k_i dx/2)/dx, as in arXiv:2512.13653 Eqs.80-81.
The midpoint current discretization is our explicit choice; the paper does not
specify enough details to claim bitwise agreement with its HILA implementation.
"""
import numpy as np


def currents(state, box, a, *, zero_policy='raise'):
    """Return J0, link Ji, radial projections and exact kinetic/gradient partitions.

    A unit projection is undefined at an exact zero field or zero link midpoint.
    'raise' is the default. Explicit 'zero' sets those projections to zero and
    records their counts and unassigned energy; it never adds an epsilon.
    """
    state=np.asarray(state,dtype=float)
    n=state.shape[0]
    if n<2 or state.shape!=(n,n,n,4) or not np.isfinite(state).all() or not np.isfinite([box,a]).all() or box<=0 or a<=0:
        raise ValueError('Require finite cubic (N,N,N,4) state, positive box and a')
    if zero_policy not in ('raise','zero'):
        raise ValueError('Unknown zero policy')
    phi=state[...,0]+1j*state[...,1]
    velocity=(state[...,2]+1j*state[...,3])/a**2
    radius=np.abs(phi);vertex_zero=radius==0
    if vertex_zero.any() and zero_policy=='raise':
        raise ValueError('Current projection undefined at exact zero vertex')
    direction=np.divide(phi,radius,out=np.zeros_like(phi),where=~vertex_zero)
    time_projection=direction.conj()*velocity
    j0=time_projection.imag;radial_time=time_projection.real
    jlink=[];radial_links=[];link_zeros=[];gradient_total=0.
    dx=box/n
    for axis in range(3):
        neighbor=np.roll(phi,-1,axis)
        midpoint=(phi+neighbor)/2
        modulus=np.abs(midpoint);zero=modulus==0
        if zero.any() and zero_policy=='raise':
            raise ValueError('Current projection undefined at exact zero link midpoint')
        direction=np.divide(midpoint,modulus,out=np.zeros_like(midpoint),where=~zero)
        derivative=(neighbor-phi)/dx
        projection=direction.conj()*derivative
        jlink.append(projection.imag);radial_links.append(projection.real)
        link_zeros.append(int(zero.sum()))
        gradient_total+=float(np.mean(np.abs(derivative)**2))/(2*a*a)
    jlink=np.array(jlink);radial_links=np.array(radial_links)
    kinetic_total=float(np.mean(np.abs(velocity)**2))/(2*a*a)
    kinetic_angular=float(np.mean(j0*j0))/(2*a*a)
    kinetic_radial=float(np.mean(radial_time*radial_time))/(2*a*a)
    gradient_angular=float(np.mean(np.sum(jlink*jlink,axis=0)))/(2*a*a)
    gradient_radial=float(np.mean(np.sum(radial_links*radial_links,axis=0)))/(2*a*a)
    partition=dict(kinetic_total=kinetic_total,kinetic_angular=kinetic_angular,
                   kinetic_radial=kinetic_radial,gradient_total=gradient_total,
                   gradient_angular=gradient_angular,gradient_radial=gradient_radial,
                   kinetic_unassigned=kinetic_total-kinetic_angular-kinetic_radial,
                   gradient_unassigned=gradient_total-gradient_angular-gradient_radial,
                   vertex_zeros=int(vertex_zero.sum()),link_midpoint_zeros=link_zeros,
                   zero_policy=zero_policy)
    return j0,jlink,partition


def mode_spectrum(j0,jlink,box,a,mass,*,dispersion='continuum'):
    """Normalized FFT mode sums are physical densities, not extensive particle counts.

    'continuum': manuscript Omega=sqrt(k²+a²m²), denominator k², K-basis.
    'lattice': stencil Omega=sqrt(K²+a²m²), denominator K², same K-basis.
    Full complex FFTs include each conjugate mode once; no extra degeneracy factor.
    The zero mode is excluded from both number sums. Its kinetic-only diagnostic
    is reported separately; spatial currents cannot reconstruct its potential.
    """
    j0=np.asarray(j0,dtype=float);jlink=np.asarray(jlink,dtype=float)
    n=j0.shape[0]
    if n<2 or j0.shape!=(n,n,n) or jlink.shape!=(3,n,n,n) or not np.isfinite(j0).all() or not np.isfinite(jlink).all():
        raise ValueError('Invalid current arrays')
    if not np.isfinite([box,a,mass]).all() or box<=0 or a<=0 or mass<0:
        raise ValueError('Invalid units or mass')
    if dispersion not in ('continuum','lattice'):
        raise ValueError('Unknown dispersion')
    dx=box/n;k1=2*np.pi*np.fft.fftfreq(n,d=dx)
    k=np.array(np.meshgrid(k1,k1,k1,indexing='ij'))
    K=2*np.sin(k*dx/2)/dx
    k2=np.sum(k*k,axis=0);K2=np.sum(K*K,axis=0)
    length=np.sqrt(K2)
    basis=np.divide(K,length[None,...],out=np.zeros_like(K),where=length[None,...]>0)
    temporal=np.fft.fftn(j0)/n**3
    spatial=np.array([np.fft.fftn(jlink[i])/n**3*np.exp(-1j*k[i]*dx/2) for i in range(3)])
    longitudinal=np.sum(basis*spatial,axis=0)
    p0=np.abs(temporal)**2;ps=np.abs(longitudinal)**2
    spatial_power=np.sum(np.abs(spatial)**2,axis=0)
    transverse=np.maximum(0.,spatial_power-ps)
    reconstruction=k2 if dispersion=='continuum' else K2
    omega=np.sqrt(reconstruction+(a*mass)**2)
    nonzero=reconstruction>0
    n0=np.divide(p0,a*omega,out=np.zeros_like(p0),where=nonzero)
    ns=np.divide(omega*ps,a*reconstruction,out=np.zeros_like(ps),where=nonzero)
    number=(n0+ns)/2
    energy=number*omega/a
    kinetic_zero=None if mass==0 else float(p0[0,0,0]/(2*a*a*mass))
    parseval0=float(abs(p0.sum()-np.mean(j0*j0)))
    parseval_spatial=float(abs(spatial_power.sum()-np.mean(np.sum(jlink*jlink,axis=0))))
    return dict(k=np.sqrt(k2),K=np.sqrt(K2),omega=omega,p0=p0,ps=ps,transverse=transverse,
                number=number,number_temporal=n0,number_spatial=ns,energy=energy,
                summary=dict(dispersion=dispersion,number_density=float(number.sum()),
                             temporal_estimator=float(n0.sum()),spatial_estimator=float(ns.sum()),
                             harmonic_axion_energy_proxy=float(energy.sum()),
                             spatial_transverse_power=float(transverse.sum()),
                             spatial_longitudinal_power=float(ps.sum()),
                             temporal_parseval_error=parseval0,spatial_parseval_error=parseval_spatial,
                             excluded_zero_mode_j0=float(j0.mean()),
                             excluded_zero_mode_kinetic_number=kinetic_zero,
                             zero_mode_policy='excluded; kinetic-only diagnostic does not recover potential'))


def shell_spectrum(modes,box):
    """Exact mode-density sums per logarithmic shell; integral equals mode total.

    Shells are i*2pi/L <= |k| < (i+1)*2pi/L, i>=1. k=0 is separate.
    Bin densities preserve actual finite-box mode multiplicities. Continuum
    k³P/(2pi²) shell-mean powers are also saved for comparison, not substituted
    for exact mode sums in sparsely populated low-k bins.
    """
    fund=2*np.pi/box;ids=np.floor(modes['k']/fund+1e-12).astype(int).ravel()
    counts=np.bincount(ids);bins=np.flatnonzero(counts);bins=bins[bins>0]
    def total(key):return np.bincount(ids,weights=modes[key].ravel(),minlength=len(counts))[bins]
    kmean=total('k')/counts[bins];dln=np.log((bins+1)/bins)
    volume=box**3
    p0mean=volume*total('p0')/counts[bins];psmean=volume*total('ps')/counts[bins]
    return dict(bin=bins,modes=counts[bins],k=kmean,k_low=bins*fund,k_high=(bins+1)*fund,
                dlnk=dln,number_shell=total('number'),number_per_dlnk=total('number')/dln,
                number_temporal_per_dlnk=total('number_temporal')/dln,
                number_spatial_per_dlnk=total('number_spatial')/dln,
                harmonic_energy_per_dlnk=total('energy')/dln,
                P0_shell_mean=p0mean,Ps_shell_mean=psmean,
                calP0=kmean**3*p0mean/(2*np.pi**2),calPs=kmean**3*psmean/(2*np.pi**2))


def normalized_yield(number_density,a,a_star,h_star,f_a=1.):
    if not np.isfinite([number_density,a,a_star,h_star,f_a]).all() or min(a,a_star,h_star,f_a)<=0 or number_density<0:
        raise ValueError('Invalid yield normalization')
    return (a/a_star)**3*number_density/(f_a*f_a*h_star)


def coherent_zero_mode(state,a,mass,*,f_a=1.,theta_limit=.3,radial_limit=.1):
    """Small-angle completion of k=0 for the finite-box toy only.

    Spatial currents cannot recover a constant angular displacement. Use the
    mean angle about the +X vacuum only after explicit harmonic-regime gates.
    This is separate from the paper-inspired nonzero-mode current estimator.
    Defect freedom must also be checked by the caller.
    """
    if min(a,mass,f_a)<=0 or not np.isfinite([a,mass,f_a]).all():
        raise ValueError('Require positive finite scale factor, mass and f_a')
    phi=state[...,0]+1j*state[...,1];radius=np.abs(phi)
    if np.any(radius==0):raise ValueError('Zero-mode phase undefined at zero field')
    theta=np.angle(phi)
    theta_max=float(np.max(abs(theta)));radial_rms=float(np.sqrt(np.mean((radius/f_a-1)**2)))
    eligible=bool(theta_max<=theta_limit and radial_rms<=radial_limit)
    j0=np.imag(np.conj(phi/radius)*(state[...,2]+1j*state[...,3])/a**2)
    rho0=float(np.mean(j0)**2/(2*a*a)+.5*mass*mass*(f_a*np.mean(theta))**2)
    return dict(eligible=eligible,theta_max=theta_max,radial_rms=radial_rms,
                number_density=rho0/mass if eligible else None,
                energy_density=rho0 if eligible else None,mean_angle=float(np.mean(theta)))

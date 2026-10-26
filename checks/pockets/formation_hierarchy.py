"""Analytic planning gates for accessible formation-field experiments.

No field evolution or grid allocation. Thin-wall/tanh quantities are
initialisation estimates, not measured tension, dynamics or survival.
"""
from dataclasses import asdict
import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn
from paper_reference import Parameters,force,hessian,terms,head_on_reflection


def false_branch(a,p):
    """Local phi minimum connected to FV at fixed instantaneous axion field.

    This is not a joint equilibrium: dV/da need not vanish. It diagnoses
    initial backreaction; imposing this branch on an oscillating field would
    be an adiabatic approximation requiring a separate timescale test.
    """
    if not np.isfinite(a):raise ValueError('Require finite axion amplitude')
    if a==0.:phi=0.
    else:
        lower=-1.
        while force(lower,a,p)[0]>0.:
            lower*=2.
            if lower < -1024:raise RuntimeError('FV branch bracket failed')
        phi=brentq(lambda x:float(force(x,a,p)[0]),lower,0.,xtol=1e-14)
    xx,xy,aa=map(float,hessian(phi,a,p))
    low=.5*(xx+aa-math.hypot(xx-aa,2*xy))
    u,va,m2,_=terms(phi,a,p)
    return dict(axion_amplitude=float(a),phi_minimum=phi,phi_force_residual=float(force(phi,a,p)[0]),
        phi_curvature=xx,axion_mass_ratio_to_unshifted=math.sqrt(float(m2))/p.mass_false,
        full_hessian_min_eigenvalue=low,adiabatic_axion_curvature=aa-xy*xy/xx,
        potential=float(u+va),interpretation='Fixed-a dilaton minimum only; not joint equilibrium or field stability')


def dirichlet_band_modes(radius,mass_false,omega_min,omega_max):
    """Reflective-sphere comparison only; finite-wall eigenmodes will differ."""
    maximum=radius*math.sqrt(omega_max**2-mass_false**2)
    x=np.linspace(1e-5,maximum,2049);rows=[]
    for ell in range(math.ceil(maximum)+1):
        y=spherical_jn(ell,x)
        intervals=np.flatnonzero(np.signbit(y[:-1])!=np.signbit(y[1:]))
        for radial,i in enumerate(intervals,1):
            root=brentq(lambda z:spherical_jn(ell,z),x[i],x[i+1])
            omega=math.sqrt((root/radius)**2+mass_false**2)
            if omega>=omega_min:
                rows.append(dict(ell=ell,radial=radial,omega=omega,angular_degeneracy=2*ell+1))
    return rows


def candidate(radius=12.,box=48.,n=192,speed=.3,omega_min=.3,omega_max=.5,dv=.05):
    if (not np.isfinite([radius,box,speed,omega_min,omega_max,dv]).all() or
        not isinstance(n,int) or n<4 or not 0<radius<box/2 or
        not 0<=speed<1 or not .1<=omega_min<=omega_max):
        raise ValueError('Invalid candidate geometry, spectrum or wall speed')
    p=Parameters(dv=dv);w=math.sqrt(8./p.lam);dx=box/n
    gamma=1/math.sqrt(1-speed*speed);sigma=math.sqrt(p.lam/2.)/6.
    mass_true=p.mass_false*math.exp(p.kappa/2.)
    edges=[head_on_reflection(omega,p.mass_false,speed) for omega in (omega_min,omega_max)]
    rest_frequencies=[v['wall_frame_frequency'] for v in edges]
    encounter_sequences=[]
    for omega in (omega_min,omega_max):
        encounters=[];incoming=omega
        for encounter in range(1,7):
            value=head_on_reflection(incoming,p.mass_false,speed)
            reflective=value['wall_frame_frequency']<mass_true
            encounters.append(dict(encounter=encounter,incoming_lab_frequency=incoming,
                **value,below_exterior_threshold=reflective))
            if not reflective:break
            incoming=value['reflected_lab_frequency']
        encounter_sequences.append(dict(initial_omega=omega,encounters=encounters,
            assumption='Constant approaching wall speed at successive head-on encounters; no dynamics or damping'))
    penetrations=[1/math.sqrt(mass_true*mass_true-f*f) if f<mass_true else None for f in rest_frequencies]
    wall_cells=2*math.atanh(.8)*w/gamma/dx
    modes=dirichlet_band_modes(radius,p.mass_false,omega_min,omega_max)
    # Conditional linear isotropic estimate, not a prescribed field pressure.
    pressure=dv+2*sigma/radius
    typical_omega=(omega_min+omega_max)/2
    rms=math.sqrt(3*pressure/(typical_omega**2-p.mass_false**2))
    # Positive reflection gap is necessary only; wall transmission, oblique
    # incidence and later spectral compression must be measured with fields.
    return dict(parameters=asdict(p),radius=radius,box=box,n=n,dx=dx,speed=speed,
        repeated_head_on_estimates=encounter_sequences,band=[omega_min,omega_max],mass_true=mass_true,wall_width_parameter=w,
        initial_tanh_10_90_cells=wall_cells,head_on_edges=edges,
        wall_frame_penetration_lengths=penetrations,
        wall_frame_penetration_cells=[v/dx if v is not None else None for v in penetrations],
        boosted_lab_penetration_cells=[v/(gamma*dx) if v is not None else None for v in penetrations],
        shortest_interior_wavelength_cells=2*math.pi/math.sqrt(omega_max**2-p.mass_false**2)/dx,
        dirichlet_fundamental_frequency=math.sqrt((math.pi/radius)**2+p.mass_false**2),
        dirichlet_band_modes=modes,dirichlet_angular_mode_count=sum(v['angular_degeneracy'] for v in modes),
        conditional_linear_isotropic_balance=dict(pressure=pressure,typical_omega=typical_omega,
            axion_rms=rms,assumptions='Linear isotropic monochromatic stress, unshifted FV mass, degenerate tension; not measured pressure or actual seed normalization'),
        degenerate_tension_reference=sigma,rc_over_radius=3*sigma/dv/radius,
        bias_over_degenerate_barrier=dv/(p.lam/64.),
        minimum_four_field_storage_gib=32.*n**3/2**30,
        initial_interface_to_image_gap=box-2*radius,
        long_lived_target_time=20*radius,
        gates=dict(initial_wall_four_cells=wall_cells>=4.,
                   initial_head_on_reflection=all(f<mass_true for f in rest_frequencies),
                   boosted_lab_penetration_four_cells=all(v is not None and v/(gamma*dx)>=4. for v in penetrations)),
        limits='Estimated initial profiles only. Boosted penetration assumes a locally steady wall and a monochromatic evanescent tail; verify actual lab gradients. No isolation, memory-capacity or survival gate is passed.')

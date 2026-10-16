"""Analytic wave, partition, projection and coordinate-normalization tests."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('pq_spectra',ROOT/'checks/pq_strings/spectra.py')
spectra=importlib.util.module_from_spec(spec);spec.loader.exec_module(spectra)


def wave(n,box,a,mass,phase,mode=(2,3,1),amplitude=1e-5,radius=1.):
    x=(np.arange(n)+.5)*box/n
    coordinates=np.meshgrid(x,x,x,indexing='ij')
    argument=sum(2*np.pi*m*c/box for m,c in zip(mode,coordinates))
    K2=sum((2*np.sin(np.pi*m/n)/(box/n))**2 for m in mode)
    omega=np.sqrt(K2+(a*mass)**2)
    theta=amplitude*np.cos(argument)*np.cos(phase)
    velocity=-amplitude*omega*np.cos(argument)*np.sin(phase)
    phi=radius*np.exp(1j*theta);p=a*a*1j*phi*velocity
    state=np.stack((phi.real,phi.imag,p.real,p.imag),axis=-1)
    return state,omega,radius**2*amplitude**2*omega/(4*a)


class PQSpectraTests(unittest.TestCase):
    def test_oblique_wave_normalization_all_phases_and_modulus_weight(self):
        for radius in (1.,2.):
            for phase in np.linspace(0,2*np.pi,9):
                state,omega,expected=wave(24,12.,1.7,.13,phase,radius=radius)
                j0,ji,_=spectra.currents(state,12.,1.7)
                modes=spectra.mode_spectrum(j0,ji,12.,1.7,.13,dispersion='lattice')
                self.assertAlmostEqual(modes['summary']['number_density']/expected,1.,delta=1e-9)
                self.assertAlmostEqual(modes['summary']['harmonic_axion_energy_proxy']/(expected*omega/1.7),1.,delta=1e-9)
                shells=spectra.shell_spectrum(modes,12.)
                self.assertAlmostEqual(np.sum(shells['number_per_dlnk']*shells['dlnk'])/expected,1.,delta=1e-9)

    def test_radial_motion_is_not_counted_as_axions_and_partition_is_exact(self):
        rng=np.random.default_rng(3);n=8;box=5.;a=1.3
        state=rng.normal(size=(n,n,n,4));state[...,:2]+=.2
        j0,ji,parts=spectra.currents(state,box,a)
        self.assertAlmostEqual(parts['kinetic_unassigned'],0.,delta=1e-12)
        self.assertAlmostEqual(parts['gradient_unassigned'],0.,delta=1e-12)
        radial=np.zeros_like(state);radial[...,0]=1.+.1*rng.random((n,n,n));radial[...,2]=rng.normal(size=(n,n,n))
        j0,ji,parts=spectra.currents(radial,box,a)
        modes=spectra.mode_spectrum(j0,ji,box,a,.1)
        self.assertEqual(modes['summary']['number_density'],0.)
        self.assertGreater(parts['kinetic_radial'],0.)
        self.assertGreater(parts['gradient_radial'],0.)

    def test_link_alignment_projection_and_parseval(self):
        n=16;box=8.;axis=(np.arange(n)+.5)*box/n
        j0=np.zeros((n,n,n));ji=np.zeros((3,n,n,n))
        ji[1]=np.cos(2*np.pi*axis/box)[:,None,None]
        modes=spectra.mode_spectrum(j0,ji,box,1.,.2)
        self.assertEqual(modes['summary']['number_density'],0.)
        self.assertAlmostEqual(modes['summary']['spatial_transverse_power'],.5,delta=1e-14)
        self.assertLess(modes['summary']['spatial_parseval_error'],1e-14)

    def test_zero_mode_is_explicit_and_undefined_projections_rejected(self):
        n=8;j0=np.ones((n,n,n))*.2;ji=np.zeros((3,n,n,n))
        modes=spectra.mode_spectrum(j0,ji,8.,2.,.1)
        self.assertEqual(modes['summary']['number_density'],0.)
        self.assertAlmostEqual(modes['summary']['excluded_zero_mode_kinetic_number'],.05)
        state=np.zeros((n,n,n,4))
        with self.assertRaises(ValueError):spectra.currents(state,8.,1.)
        _,_,parts=spectra.currents(state,8.,1.,zero_policy='zero')
        self.assertEqual(parts['vertex_zeros'],n**3)
        opposite=np.zeros_like(state);opposite[...,0]=(-1.)**np.arange(n)[:,None,None]
        with self.assertRaises(ValueError):spectra.currents(opposite,8.,1.)

    def test_coordinate_gauge_invariance_and_continuum_convergence(self):
        state,_,_=wave(16,10.,1.5,.2,.7)
        j0,ji,_=spectra.currents(state,10.,1.5)
        original=spectra.mode_spectrum(j0,ji,10.,1.5,.2,dispersion='continuum')['summary']['number_density']
        c=.5;rescaled=state.copy();rescaled[...,2:]/=c**3
        j0,ji,_=spectra.currents(rescaled,10*c,1.5/c)
        changed=spectra.mode_spectrum(j0,ji,10*c,1.5/c,.2,dispersion='continuum')['summary']['number_density']
        self.assertAlmostEqual(changed/original,1.,delta=1e-12)
        self.assertEqual(spectra.normalized_yield(original,1.5,1.,.04),spectra.normalized_yield(changed,1.5/c,1./c,.04))
        errors=[]
        for n in (16,32,64):
            state,_,expected=wave(n,12.,1.7,.13,.7)
            j0,ji,_=spectra.currents(state,12.,1.7)
            value=spectra.mode_spectrum(j0,ji,12.,1.7,.13,dispersion='continuum')['summary']['number_density']
            errors.append(abs(value/expected-1))
        self.assertGreater(errors[0]/errors[1],3.5)
        self.assertGreater(errors[1]/errors[2],3.5)


    def test_coherent_zero_mode_completion_and_harmonic_gate(self):
        a=1.7;mass=.13;amp=.1
        for phase in np.linspace(0,2*np.pi,9):
            theta=amp*np.cos(phase);velocity=-a*mass*amp*np.sin(phase)
            phi=np.exp(1j*theta);mom=a*a*1j*phi*velocity
            state=np.empty((4,4,4,4));state[:]=[phi.real,phi.imag,mom.real,mom.imag]
            zero=spectra.coherent_zero_mode(state,a,mass)
            self.assertTrue(zero['eligible'])
            self.assertAlmostEqual(zero['number_density'],.5*mass*amp**2,places=14)
        state[...,0]=np.cos(.4);state[...,1]=np.sin(.4)
        self.assertFalse(spectra.coherent_zero_mode(state,a,mass)['eligible'])
        self.assertIsNone(spectra.coherent_zero_mode(state,a,mass)['number_density'])

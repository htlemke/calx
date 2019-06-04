"""Functions for simple calculations with gaussian beams."""
import numpy as np

def waist(lam,divergence,n=1):
    """wavelength <lam> in same units as output, 
    divergence (full divergence far from focus) in radians.
    <n> is refractive index."""
    return lam/np.pi/divergence*2 /n

def wz(waist,z,lam,Msq=1.):
    """beam size at distance <z> from the waist in beam direction.
    <waist> (1/e radius in amplitude), <z> and <lam> should be same units,
    result unit will be the same."""
    M = sqrt(Msq)
    wz = w0*sqrt(M**2+M**2*(lam*z/np.pi/w0**2)**2)
    return wz

def rayleighRange(lam,waist):
    """Rayleigh range at wavelength <lam> and <waist>.
    All inputs and outputs should have same unit."""
    return np.pi*waist**2/lam

def gouyPhase(z,rayleighRange):
    """calculates the gouy phase shift close to the focus.
    distance from focus <z> and <rayleighRange> in same units, 
    result in radians"""
    return np.arctan(z/rayleighRange)



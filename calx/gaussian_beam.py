"""Functions for simple calculations with gaussian beams."""

import numpy as np


def waist(lam, divergence, n=1):
    """wavelength <lam> in same units as output,
    divergence (full divergence far from focus) in radians.
    <n> is refractive index."""
    return lam / np.pi / divergence * 2 / n


def wz(z, lam, waist=None, divergence=None, n=1, Msq=1.0):
    """beam size at distance <z> from the waist in beam direction.
    <waist> (1/e radius in amplitude), <z> and <lam> should be same units,
    result unit will be the same."""
    if waist is None:
        if divergence is None:
            raise ValueError("Either waist or divergence must be given.")
        waist = lam / np.pi / divergence * 2 / n
    M = np.sqrt(Msq)
    wz = waist * np.sqrt(M**2 + M**2 * (lam * z / np.pi / waist**2) ** 2)
    return wz


def rayleigh_range(lam, waist):
    """Rayleigh range at wavelength <lam> and <waist>.
    All inputs and outputs should have same unit."""
    return np.pi * waist**2 / lam


def gouy_phase(z, rayleighRange):
    """calculates the gouy phase shift close to the focus.
    distance from focus <z> and <rayleighRange> in same units,
    result in radians"""
    return np.arctan(z / rayleighRange)

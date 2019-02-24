import numpy as np
from scipy.constants import c,h,eV

def dThetadE(E,Q):
    """ Calculates the bragg angle derivative to energy at a certain Q and 
    photon Energy"""
    return (-*c*h*Q/eV)/(4*np.pi*E**2*np.sqrt(1-(c**2*h**2*Q**2)/(16*np.pi**2*E**2*eV**2)))

import numpy as np
from scipy import constants

def fluence(energy,area,inc_angle=np.pi/2,reflectivity=0):
    """Simple calculation of fluence at a surface of reflectivity 
    <reflectivity> with angle 
    <inc_angle> (radians) between surface plane and incoming wavevector.
    energy in J, area (normal to wavewector) in m**2.
    --> output in mJ/cm**2
    """
    return energy*1e3/area/1e4*(1-reflectivity)*np.sin(inc_angle)


def grazingIncTimeres(probeSize,probeAng, pumpAng):
  """Calculates the time resulution in a 
  grazing incidence pump probe geometry.
  This respects only purely geometric factors of the 
  surface. 
  angles in radians
  lengthscale in meters
  --> result in seconds"""
  footprint = probeSize/np.sin(probeAng)
  lenprobe = probeSize/np.tan(probeAng)
  lenpump = np.cos(pumpAng)*footprint
  return (lenpump-lenprobe)/constants.c

import numpy as np

def fluence(energy,area,inc_angle=np.pi/2,reflectivity=0):
    """Simple calculation of fluence at a surface of reflectivity 
    <reflectivity> with angle 
    <inc_angle> (radians) between surface plane and incoming wavevector.
    energy in J, area (normal to wavewector) in m**2.
    --> output in mJ/cm**2
    """
    return energy*1e3/area/1e4*(1-reflectivity)*np.sin(inc_angle)



import numpy as np
from scipy import constants


def sigr2sigx(sigr):
    return sigr / np.sqrt(2)


def sigr2fwhm(sigr):
    return sigr / np.sqrt(2) * 2.35


def dB2ratio(dB, field=False):
    """ Converts dB to a ratio, default is power ratio,
    field ratio available through keyword.
    Usages:
        field ratio: cameras signal/noise
    """
    if field:
        f = 20
    else:
        f = 10
    return 10 ** (dB / f)


def ratio2dB(ratio, field=False):
    """ Converts ratio to decibel, default is power ratio,
    field ratio available through keyword.
    Usages:
        field ratio: cameras signal/noise
    """
    if field:
        f = 20
    else:
        f = 10
    return f * np.log10(ratio)

def reccm2m(reccm):
    """ Converts wavenumer to wavelength in meters """
    return .01/reccm

def reccm2THz(reccm):
    return constants.c/0.01/1e12*reccm

def THz2reccm(THz):
    return THz/constants.c*0.01*1e12

def THz2um(THz):
    return constants.c/1e12*1e6 / THz

def um2THz(um):
    return constants.c/1e12*1e6 / um

def eV2angstrom(eV):
    return constants.c/constants.electron_volt*constants.h*1e10/eV

def angstrom2eV(angstrom):
    return constants.c/constants.electron_volt*constants.h*1e10/angstrom

import numpy as np


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

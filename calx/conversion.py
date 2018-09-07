import numpy as np

def sigr2sigx(sigr):
  return sigr/np.sqrt(2)

def sigr2fwhm(sigr):
  return sigr/np.sqrt(2)*2.35

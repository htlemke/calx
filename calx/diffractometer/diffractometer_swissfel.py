import numbers
import numpy as np
from numpy.linalg import norm
import warnings
import re

from xrayutilities import Experiment,QConversion
# package internal imports
from xrayutilities import math
# from xrayutilities import materials
# from xrayutilities import utilities
# from xrayutilities import cxrayutilities
from xrayutilities import config
from xrayutilities.exception import InputError

math.Transform
# python 2to3 compatibility
# try:
    # basestring
# except NameError:
    # basestring = str

# regular expression to check goniometer circle syntax
# directionSyntax = re.compile("[xyz][+-]")
# circleSyntaxDetector = re.compile("([xyz][+-])|(t[xyz])")
# circleSyntaxSample = re.compile("[xyzk][+-]")

class GID(Experiment):

    """
    class describing grazing incidence x-ray diffraction experiments
    the class helps with calculating the angles of Bragg reflections
    as well as it helps with analyzing measured data

    the class describes a four circle (alpha_i,azimuth,twotheta,beta)
    goniometer to help with GID experiments at the ROTATING ANODE.
    3D data can be treated with the use of linear and area detectors.
    see help self.Ang2Q

    Using this class the default sample surface orientation is determined by
    the inner most sample rotation (which is usually the azimuth motor).
    """

    def __init__(self, idir=[1,0,0], ndir=[0,0,1], **keyargs):
        """
        initialization routine for a GID Experiment class fitting a Bernina 
        diffractometer a swissfel.

        The normal convention is kept, which uses the sample z axis as the innermost rotation.
         

        idir defines the inplane reference direction (idir points into the PB
           direction at zero angles)
        ndir defines the surface normal of your sample (ndir points along the
           innermost sample rotation axis)

        Parameters
        ----------
        same as for the Experiment base class
        """
        if 'sampleor' not in keyargs:
            keyargs['sampleor'] = 'sam'

        if "qconv" not in keyargs:
            # 2S+2D goniometer
            keyargs['qconv'] = QConversion(['z-', 'x+'], ['x+', 'z-'],
                                           [0, 1, 0])

        Experiment.__init__(self, idir, ndir, **keyargs)

    def Q2Ang(self, Q, alpha_i=0, trans=True, deg=True, **kwargs):
        """
        calculate the GID angles needed in the experiment
        the inplane reference direction defines the direction were
        the reference direction is parallel to the primary beam
        (i.e. lattice planes perpendicular to the beam)

        Parameters
        ----------
         Q:          a list or numpy array of shape (3) with
                     q-space vector components

        optional keyword arguments:
         trans:      True/False apply coordinate transformation on Q
         deg:        True/Flase (default True) determines if the
                     angles are returned in radians or degrees

        Returns
        -------
        a numpy array of shape (4) with the four GID scattering angles which
        are [alpha_i, azimuth, twotheta, beta]

         alpha_i:    incidence angle to surface (at the moment always 0)
         azimuth:    sample rotation with respect to the inplane reference
                     direction
         twotheta:   scattering angle
         beta:       exit angle from surface (at the moment always 0)
        """

        for k in kwargs.keys():
            if k not in ['trans', 'deg']:
                raise Exception("unknown keyword argument given: allowed are "
                                "'trans': coordinate transformation flag, "
                                "'deg': degree-flag")

        if isinstance(Q, list):
            q = np.array(Q, dtype=np.double)
        elif isinstance(Q, np.ndarray):
            q = Q
        else:
            raise TypeError("Q vector must be a list or numpy array")

        if trans:
            q = self.Transform(q)

        if config.VERBOSITY >= config.INFO_ALL:
            print("XU.GID.Q2Ang: q = %s" % repr(q))

        if deg:
            alpha_i = np.radians(alpha_i)
            # print('converted alpha_i')

        # set parameters for the calculation
        z = self.Transform(self.ndir)  # z
        y = self.Transform(self.idir)  # y
        x = self.Transform(self.scatplane)  # x

        # check if reflection is inplane
        # if np.abs(math.VecDot(q, z)) >= 0.001:
            # raise InputError("Reflection not reachable in GID geometry (Q: %s)"
                            #  % str(q))
        #calculate perpendicular (z) and parallel components of Q 
        qz = math.VecDot(z, q)
        qx = math.VecDot(x, q)
        qy = math.VecDot(y, q)
        qxy = np.sqrt(qx**2+qy**2)

        # calculate angle to inplane reference direction, 
        # that is perpendicular to incoming beam projection on plane
        aref = np.arctan2(qx, qy)
        print(f'aref is {np.degrees(aref)} degrees')

        # calculate scattering angle
        qa = math.VecNorm(q)
        tth = 2. * np.arcsin(qa / 2. / self.k0)
        # calculate "phi" rotation that brings into ewald sphere
        a_ewald = np.arcsin(
            qa**2/self.k0/2/np.cos(alpha_i)/qxy + qz/qxy*np.tan(alpha_i)
            )
        # print(f'a_ewald is {np.degrees(a_ewald)} degrees')
        # print(f'alpha_i is {np.degrees(alpha_i)} degrees')

        azimuth = np.pi / 2 + aref + a_ewald
        # print(f'az is {np.degrees(azimuth)} degrees')
        q_lab = self.Ang2Q.transformSample2Lab(q,np.degrees(alpha_i),np.degrees(azimuth))
        k_out = q_lab + self.k0*self.idir
        # print('k out comparison',norm(self.k0*self.idir+q_lab),self.k0)
        e_out = q_lab/self.k0 + self.idir
        # print(f'ql is {q_lab} q is {q}')
        assert len(self.Ang2Q.detectorAxis) == 2, "we need 2 detector axes for this"
        ax0 = math.getVector(self.Ang2Q.detectorAxis[0])
        ax1 = math.getVector(self.Ang2Q.detectorAxis[1])
        a1 = np.pi/2 - np.arccos(math.VecDot(ax0,e_out))
        # print(ax0,ax1,e_out)
        ax1r = math.VecCross(e_out,ax0)
        e_out0 = math.rotarb(e_out,ax1r,-a1,deg=False)
        a0 = - np.arccos(math.VecDot(e_out0,self.idir))



        if deg:
            ang = [ np.degrees(azimuth), np.degrees(a0), np.degrees(a1)]
        else:
            ang = [azimuth, a0, a1]

        if config.VERBOSITY >= config.INFO_ALL:
            print("XU.GID.Q2Ang: [ai,azimuth,tth,beta] = %s \n difference to "
                  "inplane reference which is %5.2f" % (str(ang), aref))

        return ang

    def Ang2Q(self, ai, phi, tt, beta, **kwargs):
        """
        angular to momentum space conversion for a point detector. Also see
        help GID.Ang2Q for procedures which treat line and area detectors

        Parameters
        ----------
         ai,phi,tt,beta: sample and detector angles as numpy array, lists or
                         Scalars must be given. All arguments must have the
                         same shape or length. However, if one angle is always
                         the same its enough to give one scalar value.

        **kwargs:   optional keyword arguments
            delta:  giving delta angles to correct the given ones for
                    misalignment delta must be an numpy array or list of
                    length 4. Used angles are than ai,phi,tt,beta - delta
            UB:     matrix for conversion from (hkl) coordinates to Q of sample
                    used to determine not Q but (hkl)
                    (default: identity matrix)
            wl:     x-ray wavelength in angstroem (default: self._wl)
            deg:    flag to tell if angles are passed as degree (default: True)

        Returns
        -------
        reciprocal space positions as numpy.ndarray with shape ( * , 3 )
        where * corresponds to the number of points given in the input
        """
        # dummy function to have some documentation string available
        # the real function is generated dynamically in the __init__ routine
        pass
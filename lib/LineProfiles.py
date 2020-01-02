import numpy as np
from numba import njit


#Define the constant
cSqrtPI = np.sqrt(np.pi)
cSqrtLN2 = np.sqrt(np.log(2))
FinalConst = cSqrtPI*cSqrtLN2

cLn2 = 0.6931471805599
cSqrtLn2divSqrtPi = 0.469718639319144059835

A = np.array([-1.2150, -1.3509, -1.2150, -1.3509])
B = np.array([1.2359, 0.3786, -1.2359, -0.3786])
C = np.array([-0.3085, 0.5906, -0.3085, 0.5906])
D = np.array([0.0210, -1.1858, -0.0210, 1.1858])


def PROFILE_VOIGT(sg0,GamD,Gam0,sg):
    """
    # Voigt profile based on astropy module.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   Gam0: Speed-averaged line-width in cm-1 (Input).
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    X = (sg- sg0)*2*cSqrtLN2/GamD
    X = np.atleast_1d(X)[..., np.newaxis]
    Y = Gam0*cSqrtLN2/GamD
    Y = np.atleast_1d(Y)[..., np.newaxis]
    V = np.sum((C*(Y - A) + D*(X-B))/((Y-A)**2 + (X-B)**2), axis=-1)
    return (Gam0*FinalConst/GamD)*V


def PROFILE_LORENTZ(sg0,GamD,Gam0,sg):
    """
    # Lorentz profile.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   Gam0: Speed-averaged line-width in cm-1 (Input).
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return Gam0/(pi*(Gam0**2+(sg-sg0)**2))

def PROFILE_DOPPLER(sg0,GamD,Gam0,sg):
    """
    # Doppler profile.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   Gam0: Speed-averaged line-width in cm-1 (Input).
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    return cSqrtLn2divSqrtPi*np.exp(-cLn2*((sg-sg0)/GamD)**2)/GamD

import numpy as np
from numba import njit
from scipy.special import wofz


#Define the constant
cSqrtPI = np.sqrt(np.pi)
cSqrtLN2 = np.sqrt(np.log(2))
FinalConst = cSqrtPI*cSqrtLN2

cLn2 = 0.6931471805599
cSqrtLn2divSqrtPi = 0.469718639319144059835


def PROFILE_VOIGT(sg0,GamD,Gam0,sg):
    """
    # Voigt profile based on astropy module.
    # Input parameters:
    #   sg0: Unperturbed line position in cm-1 (Input).
    #   GamD: Doppler HWHM in cm-1 (Input)
    #   Gam0: Speed-averaged line-width in cm-1 (Input).
    #   sg: Current WaveNumber of the Computation in cm-1 (Input).
    """
    sg = sg - sg0
    sigma = GamD/np.sqrt(2 * np.log(2))
    return np.real(wofz(( sg + 1j*Gam0)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)





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

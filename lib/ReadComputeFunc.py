import glob
import numpy as np
from .Constants import *
from bisect import bisect
from numba import njit, prange

def ReadData(MoleculeName, Location="data/"):
    """
    Parameters
    ----------------------------------------------------------------------------
    MoleculeName: The name of the molecule for which the cross-section is being
                  generated.

    Location: Location where the datis looked at
    """
    if len(Location)>0.5:
        if Location[-1]=="/":
            DataFiles = glob.glob(Location+"*.data")
        else:
            DataFiles = glob.glob(Location+"/*.data")
    else:
        DataFiles = glob.glob("data/*.data")
    DataFiles = np.array(DataFiles)
    AllMolNames = np.array([Item.split("/")[-1][:-5].upper() for Item in DataFiles])
    SelectIndex = AllMolNames == MoleculeName.upper()

    assert np.sum(SelectIndex)==1, "No files found for the molecule"

    DataContent = open(DataFiles[SelectIndex][0],'r').readlines()

    return DataContent

#@njit(parallel=True, nopython=False, nogil=True)
def GenerateCrossSection(Omegas, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, LINE_PROFILE, Params):

    P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor = Params
    i = 0
    Xsect = np.zeros(len(Omegas))
    NLINES = len(LineCenterDB)
    for i in prange(NLINES):
        LineCenter = LineCenterDB[i]
        LowerStateEnergy = LowerStateEnergyDB[i]


        #Calculate the line shift

        #Calculate the temperature dependence of the gamma
        Gamma0 = GammaSelf[i]*P/Pref*(Tref/Temp)**TempRatioPower[i]

        #The Doppler Broadening coefficients
        GammaD = np.sqrt(2*cBolts*Temp*np.log(2)/m/cc**2)*LineCenterDB[i]

        #Scale the line intensity with the temperature
        ch = np.exp(-const_R*LowerStateEnergy/Temp)*(1-np.exp(-const_R*LineCenter/Temp))
        zn = np.exp(-const_R*LowerStateEnergy/Tref)*(1-np.exp(-const_R*LineCenter/Tref))
        LineIntensity = LineIntensityDB[i]*SigmaTref/SigmaT*ch/zn

        OmegaWingF = max(OmegaWing,OmegaWingHW*Gamma0,OmegaWingHW*GammaD)

        #Find the range to calculate the value:::
        BoundIndexLower = bisect(Omegas,LineCenterDB[i]-OmegaWingF)
        BoundIndexUpper = bisect(Omegas,LineCenterDB[i]+OmegaWingF)

        lineshape_vals = LINE_PROFILE(LineCenter,GammaD,Gamma0,Omegas[BoundIndexLower:BoundIndexUpper])
        Xsect[BoundIndexLower:BoundIndexUpper] += factor*LineIntensity*lineshape_vals
        i+=1
    return Xsect

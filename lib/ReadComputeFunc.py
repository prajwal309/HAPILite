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

    Data = open(DataFiles[SelectIndex][0],'r').readlines()
    MoleculeNumberDB = int(Data[0][0:2].replace(" ",""))
    IsoNumberDB = int(Data[0][2:3].replace(" ",""))

    #Read the important parameters from the file
    LineCenterDB = np.array([float(Item[3:15]) for Item in Data])
    LineIntensityDB = np.array([float(Item[16:26]) for Item in Data])
    LowerStateEnergyDB = np.array([float(Item[46:56].replace("-","")) for Item in Data])
    GammaSelf = np.array([float(Item[40:46]) for Item in Data])

    #Temperature Dependence of Gamma0
    TempRatioPower = np.array([float(Item[55:59]) for Item in Data])

    #The value of the error values. Keep them as string to preserve 0 at the beginning.
    ErrorArray = np.array([Item[127:133] for Item in Data])

    return MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, ErrorArray



def GenerateCrossSection(Omegas, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, LINE_PROFILE, Params):
    """
    Parameters
    --------------------------------------
    Omega: The wavegrid of the WaveNumber
    LineCenterDB: The line center from the HITRAN database
    LineIntensityDB: The line intensity from the HITRAN database
    LowerStateEnergyDB: LowerStateEnergyDB from lower state energy
    GammaSelf: The self broadening parameters
    TempRatioPower: Parameter for scaling the self broadening parameters in temperature
    LINE_PROFILE: PROFILE_VOIGT, PROFILE_DOPPLER, PROFILE_LORENTZ are the
    """
    P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor = Params
    i = 0
    Xsect = np.zeros(len(Omegas))
    NLINES = len(LineCenterDB)
    for i in prange(NLINES):
        LineCenter = LineCenterDB[i]
        LowerStateEnergy = LowerStateEnergyDB[i]

        #Calculate the line shift for hydrogen and helium

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
        BoundIndexLower = bisect(Omegas,LineCenter-OmegaWingF)
        BoundIndexUpper = bisect(Omegas,LineCenter+OmegaWingF)

        lineshape_vals = LINE_PROFILE(LineCenter,GammaD,Gamma0,Omegas[BoundIndexLower:BoundIndexUpper])
        Xsect[BoundIndexLower:BoundIndexUpper] += factor*LineIntensity*lineshape_vals
        i+=1
    return Xsect

#import standard libraries
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from bisect import bisect

#import custom functions
from lib.Constants import *
from lib.LineProfiles import PROFILE_DOPPLER, PROFILE_LORENTZ, PROFILE_VOIGT
from lib.PartitionFunction import BD_TIPS_2017_PYTHON
from lib.MolecularMass import GetMolecularMass
from lib.ReadComputeFunc import ReadData, GenerateCrossSection

from numba import jit

def CalcCrossSection(MoleculeName, DataPath = "", Temp=296, P=1, WN_Grid=np.arange(0,15000,0.01), OmegaWing=25., OmegaWingHW=75.0, Profile="VOIGT", NCORES=1):
    """
    Parameters
    ----------------------------------------------------------------------------
    MoleculeNameName: Name of the molecule for which the cross section is to be calculated.

    DataLocation: location where the HITRAN data is located.

    Temp: The temperature at which the the cross section is to be calculated

    P: The pressure at which the cross-section is to be calculated

    WN_Grid: The range of wavenumber over which to calculcate the cross-section

    Profile: Can be Doppler, Lorentz or Voigt.
    ----------------------------------------------------------------------------
    """

    #Read the data from the provided path
    if len(DataPath)>0.5:
        Data = ReadData(MoleculeName, DataPath)
    else:
        Data = ReadData(MoleculeName)   #By default look at the data folder


    MoleculeNumberDB = int(Data[0][0:2].replace(" ",""))
    IsoNumberDB = int(Data[0][2:3].replace(" ",""))

    #Read the important parameters from the file
    LineCenterDB = np.array([float(Item[3:15]) for Item in Data])
    LineIntensityDB = np.array([float(Item[16:26]) for Item in Data])
    LowerStateEnergyDB = np.array([float(Item[46:56]) for Item in Data])
    GammaSelf = np.array([float(Item[40:46]) for Item in Data])

    #ErrorIndex = np.array([int(Item[]) for Item in Data ])

    #Temperature Dependence of Gamma0
    TempRatioPower = np.array([float(Item[56:59]) for Item in Data])

    NLINES = len(LineCenterDB)
    #Units --- Not HITRAN Units
    factor = (P/9.869233e-7)/(cBolts*Temp) # CGS np.arange(0,10000,0.01)

    #Calculate the partition Function
    SigmaT = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Temp)
    SigmaTref = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Tref)

    if "DOPPLER" in Profile.upper():
        print("Using Doppler Profile")
        LINE_PROFILE = PROFILE_DOPPLER
    elif "LORENTZ" in Profile.upper():
        print("Using Lorentz Profile")
        LINE_PROFILE = PROFILE_LORENTZ
    elif "VOIGT" in Profile.upper():
        print("Using the Voigt Profile")
        LINE_PROFILE = PROFILE_VOIGT
    m = GetMolecularMass(MoleculeNumberDB,IsoNumberDB)*cMassMol*1000.0

    #Originate in the line cross-section
    Omegas = WN_Grid[:]


    Tolerance = 3.0
    SelectIndex = np.logical_and(LineCenterDB>min(Omegas)-Tolerance, LineCenterDB<max(Omegas)+Tolerance)



    #Now slice the range
    LineCenterDB = LineCenterDB[SelectIndex]
    LineIntensityDB = LineIntensityDB[SelectIndex]
    LowerStateEnergyDB = LowerStateEnergyDB[SelectIndex]
    GammaSelf = GammaSelf[SelectIndex]
    TempRatioPower = TempRatioPower[SelectIndex]


    #This can be made faster by multiple threading

    print(NCORES)
    if NCORES==1 or NLINES<5000:
        Params = [P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor]
        Xsect = GenerateCrossSection(Omegas, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, LINE_PROFILE, Params)
        return Xsect
    else:
        if NCORES == -1:
            NUMCORES = mp.cpu_count()
        elif NCORES>1.99:
            NUMCORES = int(NCORES)

        CPU_Pool = mp.Pool(NUMCORES)
        Tasks = []

        Params = [P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor]
        for i in range(NUMCORES):
            StartIndex = int(i*NLINES/NUMCORES)
            if i==NUMCORES-1:
                StopIndex = -1
            else:
                StopIndex = int((i+1)*NLINES/NUMCORES)

            #Slicing the data
            LineCenterDB_Slice = LineCenterDB[StartIndex:StopIndex]
            LineIntensityDB_Slice = LineIntensityDB[StartIndex:StopIndex]
            LowerStateEnergyDB_Slice = LowerStateEnergyDB[StartIndex:StopIndex]
            GammaSelf_Slice = GammaSelf[StartIndex:StopIndex]
            TempRatioPower_Slice = TempRatioPower[StartIndex:StopIndex]

            Tasks.append(CPU_Pool.apply_async(GenerateCrossSection, (Omegas, LineCenterDB_Slice, LineIntensityDB_Slice,
                                              LowerStateEnergyDB_Slice, GammaSelf_Slice, TempRatioPower_Slice, LINE_PROFILE, Params)))
        CPU_Pool.close()
        CPU_Pool.join()

        Xsect = np.zeros(len(Omegas))
        for task in Tasks:
            Xsect+=task.get()
        return Xsect




    i = 0
    Xsect = np.zeros(len(Omegas))
    NLINES = len(LineCenterDB)
    while i<NLINES:
        LineCenter = LineCenterDB[i]
        LowerStateEnergy = LowerStateEnergyDB[i]


        #Calculate the temperature dependence of the gamma
        Gamma0 = GammaSelf[i]
        Value = Gamma0*P/Pref*(Tref/Temp)**TempRatioPower[i]
        Gamma0=Value


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





WaveNumber = np.arange(0,10000,0.001)

import time

StartTime = time.time()
CrossSection =  CalcCrossSection("CO2",Temp=1000.0,WN_Grid=WaveNumber, Profile="Doppler", NCORES=-1)
StopTime = time.time()

print("The difference between the start and the stop time is given by:", StopTime - StartTime)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(WaveNumber, CrossSection, "k-")
plt.title("L-HAPI")
plt.show()

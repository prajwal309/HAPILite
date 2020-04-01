#import standard libraries
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from bisect import bisect

#import custom functions
from lib.Constants import *
from lib.LineProfiles import PROFILE_DOPPLER, PROFILE_LORENTZ, PROFILE_VOIGT, PROFILE_PSEUDOVOIGT
from lib.PartitionFunction import BD_TIPS_2017_PYTHON
from lib.MolecularMass import GetMolecularMass
from lib.ReadComputeFunc import GenerateCrossSection
from lib.ErrorParser import MapError

def CalcCrossSection(Database, DataPath = "", Temp=296, P=1, Broadening="Self", \
    WN_Grid=np.arange(333.33,33333.33,0.01), OmegaWing=0.0, OmegaWingHW=75.0, \
    Profile="VOIGT", NCORES=-1):
    """
    Parameters
    ----------------------------------------------------------------------------
    MoleculeNameName: Name of the molecule for which the cross section is to be calculated.

    DataLocation: location where the HITRAN data is located.

    Temp: The temperature at which the the cross section is to be calculated

    P: The pressure at which the cross-section is to be calculated

    Broadening: Broadening can either be self or air

    OmegaWing: The width of the WaveNumber for calculating the cross-section

    OmegaWingHW: The full width half maximum for calculating the cross-section

    WN_Grid: The range of wavenumber over which to calculcate the cross-section. The default value is roughly the range of the JWST spectral range.

    Profile: Can be Doppler, Lorentz or Voigt.

    NCORES: Use the number of cores for calculating the cross-section.
    ----------------------------------------------------------------------------
    """


    MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, \
    GammaSelf, GammaAir, DeltaAir, TempRatioPower, ErrorArray = Database

    #Units --- HITRAN units should be true
    factor = 1.0  #(P/9.869233e-7)/(cBolts*Temp) #

    #Calculate the partition Function
    SigmaT = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Temp)
    SigmaTref = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Tref)

    if "DOPPLER" in Profile.upper():
        #print("Using Doppler Profile")
        LINE_PROFILE = PROFILE_DOPPLER
    elif "LORENTZ" in Profile.upper():
        #print("Using Lorentz Profile")
        LINE_PROFILE = PROFILE_LORENTZ
    elif "PSEUDOVOIGT" in Profile.upper():
        #print("Using the Pseudo-Voigt Profile")
        LINE_PROFILE = PROFILE_PSEUDOVOIGT
    elif "VOIGT" in Profile.upper():
        #print("Using the Voigt Profile")
        LINE_PROFILE = PROFILE_VOIGT
    else:
        raise("No such profile found.")

    m = GetMolecularMass(MoleculeNumberDB,IsoNumberDB)*cMassMol*1000.0

    #Originate in the line cross-section
    Omegas = WN_Grid[:]

    Tolerance = 25.0
    SelectIndex = np.logical_and(LineCenterDB>min(Omegas)-Tolerance, LineCenterDB<max(Omegas)+Tolerance)



    #Now slice the range
    LineCenterDB = LineCenterDB[SelectIndex]
    LineIntensityDB = LineIntensityDB[SelectIndex]
    LowerStateEnergyDB = LowerStateEnergyDB[SelectIndex]
    GammaSelf = GammaSelf[SelectIndex]
    GammaAir =  GammaAir[SelectIndex]
    DeltaAir = DeltaAir[SelectIndex]
    TempRatioPower = TempRatioPower[SelectIndex]

    #check for nans in the values
    NLINES = len(LineCenterDB)
    Params = [P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor, Broadening]

    #how may threads to implement
    if NLINES<100 or NCORES==1:
        print("Using single core of generating cross-section")
        Xsect = GenerateCrossSection(Omegas, LineCenterDB, LineIntensityDB, LowerStateEnergyDB,\
                                     GammaSelf, GammaAir, DeltaAir, TempRatioPower, LINE_PROFILE,\
                                     Params)
        return Xsect
    else:
        if NCORES == -1:
            NUMCORES = mp.cpu_count()
        elif NCORES>1.99:
            NUMCORES = int(NCORES)

        #print("Using %d cores for generating cross-section" %NUMCORES)
        CPU_Pool = mp.Pool(NUMCORES)
        Tasks = []
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
            GammaAir_Slice =  GammaAir[StartIndex:StopIndex]
            DeltaAir_Slice = DeltaAir[StartIndex:StopIndex]
            TempRatioPower_Slice = TempRatioPower[StartIndex:StopIndex]

            Tasks.append(CPU_Pool.apply_async(GenerateCrossSection, (Omegas, LineCenterDB_Slice, \
                                                LineIntensityDB_Slice, LowerStateEnergyDB_Slice,\
                                                GammaSelf_Slice, GammaAir_Slice, DeltaAir_Slice,\
                                                TempRatioPower_Slice, LINE_PROFILE, Params)))

        CPU_Pool.close()
        CPU_Pool.join()

        Xsect = np.zeros(len(Omegas))
        for task in Tasks:
            Xsect+=task.get()
        return Xsect
    return Xsect




def CalcCrossSectionWithError(Database, DataPath = "", Temp=296.0, P=1.0, Broadening="Self", \
                              WN_Grid=np.arange(0,15000,0.01), OmegaWing=0.0, OmegaWingHW=75.0,\
                              Profile="VOIGT", NCORES=-1, Err="0Sig"):
        """
        Parameters
        ----------------------------------------------------------------------------
        MoleculeNameName: Name of the molecule for which the cross section is to be calculated.

        DataLocation: location where the HITRAN data is located.

        Temp: The temperature at which the the cross section is to be calculated

        P: The pressure at which the cross-section is to be calculated

        Broadening:The pressure broadening can be self or other method.

        OmegaWing: The width of the WaveNumber for calculating the cross-section

        OmegaWingHW: The full width half maximum for calculating the cross-section

        WN_Grid: The range of wavenumber over which to calculcate the cross-section

        Profile: Can be Doppler, Lorentz or Voigt.

        NCORES: Use the number of cores for calculating the cross-section.

        Err: The error can be 0Sig, +1Sig, +2Sig, -1Sig, -2Sig
        ----------------------------------------------------------------------------
        """

        #Considering the error
        if Err.upper() == "0SIG":
            ErrorSTD = 0.0
        elif "-1SIG" in Err.upper():
            ErrorSTD = -1.0
        elif "-2SIG" in Err.upper():
            ErrorSTD = -2.0
        elif "1SIG" in Err.upper():
            ErrorSTD = 1.0
        elif "2SIG" in Err.upper():
            ErrorSTD = 2.0
        else:
            print("The value of Err is given by:", Err)
            raise Exception("Error in HAPILite.py. The values allowed for Err are 0SIG, 1SIG, 2SIG, -1SIG, -2SIG.")

        #Unpack the database
        MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, \
        GammaSelf, GammaAir, DeltaAir, TempRatioPower, ErrorArray = Database
        ErrorValues = MapError(ErrorArray)

        #Replace nan with zero
        NLINES = len(LineCenterDB)
        #Units --- Not HITRAN Units
        factor = 1.0#(P/9.869233e-7)/(cBolts*Temp) #

        #Calculate the partition Function
        SigmaT = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Temp)
        SigmaTref = BD_TIPS_2017_PYTHON(MoleculeNumberDB,IsoNumberDB,Tref)


        if "DOPPLER" in Profile.upper():
            #print("Using Doppler Profile")
            LINE_PROFILE = PROFILE_DOPPLER
        elif "LORENTZ" in Profile.upper():
            #print("Using Lorentz Profile")
            LINE_PROFILE = PROFILE_LORENTZ
        elif "PSEUDOVOIGT" in Profile.upper():
            #print("Using the Pseudo-Voigt Profile")
            LINE_PROFILE = PROFILE_PSEUDOVOIGT
        elif "VOIGT" in Profile.upper():
            #print("Using the Voigt Profile")
            LINE_PROFILE = PROFILE_VOIGT
        else:
            raise("No such profile found.")


        m = GetMolecularMass(MoleculeNumberDB,IsoNumberDB)*cMassMol*1000.0

        #Originate in the line cross-section
        Omegas = WN_Grid[:]


        Tolerance = max(OmegaWing, OmegaWingHW*0.25)
        SelectIndex = np.logical_and(LineCenterDB>min(Omegas)-Tolerance, LineCenterDB<max(Omegas)+Tolerance)


        #Now slice the range
        LineCenterDB = LineCenterDB[SelectIndex]+ErrorValues[SelectIndex,0]*ErrorSTD
        LineIntensityDB = LineIntensityDB[SelectIndex]*(1.0+ErrorValues[SelectIndex,1]*ErrorSTD)
        LowerStateEnergyDB = LowerStateEnergyDB[SelectIndex]
        GammaSelf = GammaSelf[SelectIndex]*(1.0+ErrorValues[SelectIndex,3]*ErrorSTD)
        GammaAir =  GammaAir[SelectIndex]*(1.0+ErrorValues[SelectIndex,2]*ErrorSTD)
        DeltaAir = DeltaAir[SelectIndex]*(1.0+ErrorValues[SelectIndex,5]*ErrorSTD)
        TempRatioPower = TempRatioPower[SelectIndex]*(1.0+ErrorValues[SelectIndex,4]*ErrorSTD)

        Params = [P, Temp, OmegaWing, OmegaWingHW, m, SigmaT, SigmaTref, factor, Broadening]
        #how may threads to implement
        if NCORES==1:
            #print("Using single core of generating cross-section")
            Xsect = GenerateCrossSection(Omegas, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, LINE_PROFILE, Params)
            return Xsect
        else:
            if NCORES == -1:
                NUMCORES = mp.cpu_count()
            elif NCORES>1.99:
                NUMCORES = int(NCORES)

            #print("Using %d cores for generating cross-section" %NUMCORES)
            CPU_Pool = mp.Pool(NUMCORES)
            Tasks = []
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
                GammaAir_Slice =  GammaAir[StartIndex:StopIndex]
                DeltaAir_Slice = DeltaAir[StartIndex:StopIndex]
                TempRatioPower_Slice = TempRatioPower[StartIndex:StopIndex]

                Tasks.append(CPU_Pool.apply_async(GenerateCrossSection, (Omegas, LineCenterDB_Slice, \
                            LineIntensityDB_Slice, LowerStateEnergyDB_Slice,\
                            GammaSelf_Slice, GammaAir_Slice, DeltaAir_Slice,\
                            TempRatioPower_Slice, LINE_PROFILE, Params)))

            CPU_Pool.close()
            CPU_Pool.join()

            Xsect = np.zeros(len(Omegas), dtype=np.float32)
            for task in Tasks:
                Xsect+=task.get()
            return Xsect

        return Xsect.astype(np.float32)

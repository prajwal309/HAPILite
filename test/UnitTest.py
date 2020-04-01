#author: Prajwal Niraula
#Description: Module to assert that the output from HAPILite is the same as

import numpy as np
import time
import matplotlib.pyplot as plt

from HAPILite import CalcCrossSection
from lib.ReadComputeFunc import ReadData


def SpeedTest():
    #Test how long generating the cross-section takes.
    pass

def ProfileTest():

    "This function is supposed to test  voigt profile."

    Sg = np.linspace(990,1010, 0.01)
    #Use a range of values
    PseudoVoigtProfile = lib.LineProfile.PROFILE_PSEUDOVOIGT(Sg,0.2,0.2,1000)

    plt.figure()
    plt.plot(Sg, PseudoVoigtProfile, "ko")
    plt.show()




def VoigtTest_Air():
    Molecule = "CO"
    TempValue = 1000.0
    P_Value = 0.1#0.00001


    OmegaWingValue = 500.0
    OmegaRangeValue = [2086.32193,2093.410992]
    WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.01, 0.01)
    Env= {'T': TempValue, 'p': P_Value}


    hapi.db_begin('TestData')

    #Load the data
    nu_Hapi, abs_Hapi = np.loadtxt("CO_hapi_data.npy")


    CrossSectionDoppler =  CalcCrossSection(Database,Temp=TempValue,P = P_Value, Broadening="Self", WN_Grid=WaveNumber, Profile="Doppler", OmegaWing=OmegaWingValue, OmegaWingHW=0.0, NCORES=1)
    CrossSectionVoigt =  CalcCrossSection(Database,Temp=TempValue, P = P_Value, Broadening="Air", WN_Grid=WaveNumber, Profile="PseudoVoigt", OmegaWing=OmegaWingValue, OmegaWingHW=0.0, NCORES=1)
    TimeTakenHAPILite = time.time() - StartTime



if __main__ == "__user__":
    ProfileTest()
    #Check if the values of the light curves are

import glob
import numpy as np
from CrossSectionFunctions import GetWaveNumbers, SymplecticInterpolation
import matplotlib.pyplot as plt
import os
#from numba import jit, njit


#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])    #Step size of the temperature
TempStop = float(Values[1])     #Step size of the temperature
TempStep = float(Values[2])     #Step size of the temperature

expP_Start = float(Values[3])
expP_Stop = float(Values[4])
expP_Step = float(Values[5])                                                       #Step size of the pressure


Broadener = Values[6].replace(" ","")                                           #Broadening either self or air at this point
OmegaWidth = float(Values[7])                                                   #Consider the omegawidth -- how far the lines have to be considered
LowWavelength = float(Values[8])                                                #Shortest Wavelength coverage range
HighWavelength = float(Values[9])                                               #Longest Wavelength coverage range
WN_Resolution = float(Values[10])                                        #Resolution of the Wave Number
LineShapeProfile = Values[11].replace(" ","")                                    #Voigt profile by default
NumChunks = int(Values[12].replace(" ",""))                                      #Number of chunks for the wavenumber
MoleculeList = Values[13].split(",")                                             #Get the list of Molecular species
Cores = int(Values[14])                                                          #Number of cores to be used

MoleculeList = [Item.replace(" ", "").replace("\t","") for Item in MoleculeList]


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)                    #Temperature in K
expP_Range = np.arange(expP_Start, expP_Stop-expP_Step, -expP_Step)             #Pressure in log(P) atm

WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)            #in per cm
WaveNumberRange = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)


print("The range of temperature is given by::", TempRange)
print("The range of pressure is given by::", expP_Range)

input("Wait here....")


WaveNumberStart = 1./(HighWavelength*1.e-7)         #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)            #in per cm
WaveNumberRange = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)




#Now get the assign the resolution values
Resolution = 25000

#Low resolution wavelength
Wavelength_LR, WaveNumber_LR = GetWaveNumbers(LowWavelength, HighWavelength, Resolution)

print("The range of the wavelength is given by::", Wavelength_LR)
Folder2Save = "R"+str(Resolution)

if not(os.path.exists(Folder2Save)):
    os.system("mkdir %s" %(Folder2Save))

np.savetxt(Folder2Save+"/Temperature.txt", TempRange)
np.savetxt(Folder2Save+"/exp_Pressure.txt", expP_Range)
np.savetxt(Folder2Save+"/WaveLength.txt", Wavelength_LR)


#converting to nm
Wavelength_LR*=1e7

plt.figure()
plt.plot(Wavelength_LR,Wavelength_LR, "ko")
plt.show()


LengthWaveNumber = len(np.arange(WaveNumberStart, WaveNumberStop+WN_Resolution, WN_Resolution))
WaveNumberRanges = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)


HR_WaveNumber = WaveNumberRanges
LR_WaveNumber = np.linspace(min(HR_WaveNumber), max(HR_WaveNumber),Resolution)


for Molecule in MoleculeList:

    #Read the molecule name
    Location = "DataMatrix/"+Molecule+".npy"
    print("The Location is given by::", Location)
    SigmaMatrix = np.load(Location)
    print("The shape of the cross-section is given by::", np.shape(SigmaMatrix))
    input("Wait here....")
    #which molecule to consider


    #Initiate a molecule matrix
    MoleculeMatrix = np.zeros((len(TempRange),len(expP_Range),Resolution))

    for TempCounter in range(len(TempCount)):
        for PCounter in range((len(expP_Range))):
            print(TempCounter, PCounter)
    #Convert WaveNumbers to Wavelengths
    Selected_HR_Wavelength = 1./Selected_WN*1e7
    Selected_LR_Wavelength = 1./SelectedWaveNumber*1e7

    MinWavelength  = min(Selected_HR_Wavelength)
    MaxWavelength = max(Selected_HR_Wavelength)

    #Reverse the sigma when assigning the value
    AssignIndexStart = np.argmin(np.abs(MinWavelength-Wavelength_LR))
    AssignIndexStop = AssignIndexStart + len(InterpolatedValues)


    MoleculeMatrix[OuterCounter, InnerCounter, AssignIndexStart:AssignIndexStop] = InterpolatedValues[::-1]

    if 1==2:        #Printing in flagging
        plt.figure(figsize=(18,8),dpi=300)
        plt.plot(Selected_HR_Wavelength, GenCrossSection, "k.-", label="")
        plt.plot(Selected_LR_Wavelength, InterpolatedValues, "r+:", label="Interpolation")
        SaveName = "Figures/Case_"+str(OuterCounter*1000+InnerCounter)+".png"
        TitleText = "Temp_"+str(int(TemperatureValue)) + "_Pressure_"+str((PressureValue))
        print("The Text Value is given by::", TitleText)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Cross-Section")
        plt.title(TitleText)
        plt.tight_layout()
        plt.savefig("Figures/"+TitleText+".png")
        plt.show()
        plt.close('all')

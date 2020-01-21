import glob
import numpy as np
from CrossSectionFunctions import GetWaveNumbers, SymplecticInterpolation
import matplotlib.pyplot as plt
import os

#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])    #Step size of the temperature
TempStop = float(Values[1])     #Step size of the temperature
TempStep = float(Values[2])     #Step size of the temperature

expP_Start = float(Values[3])   #The largest log10(pressure) in atm
expP_Stop = float(Values[4])    #The smallest log10(pressure) in atm
expP_Step = float(Values[5])    #Step size of the pressure


Broadener = Values[6].replace(" ","")                                           #Broadening either self or air at this point
OmegaWidth = float(Values[7])                                                   #Consider the omegawidth -- how far the lines have to be considered
LowWavelength = float(Values[8])                                                #Shortest Wavelength coverage range
HighWavelength = float(Values[9])                                               #Longest Wavelength coverage range
WN_Resolution = float(Values[10])                                               #Resolution of the Wave Number
LineShapeProfile = Values[11].replace(" ","")                                   #Voigt profile by default
MoleculeList = Values[12].split(",")                                             #Get the list of Molecular species
Cores = int(Values[13])
Error = Values[14].replace("\t","")

MoleculeList = np.array([Item.replace(" ", "").replace("\t","") for Item in MoleculeList])

print(MoleculeList)
input("Wait here...")


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)                    #Temperature in K
expP_Range = np.arange(expP_Start, expP_Stop-expP_Step, -expP_Step)             #Pressure in log(P) atm


#Define the wavenumber range values...
WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)               #in per cm
WaveNumberRange = np.arange(WaveNumberStart, WaveNumberStop, WN_Resolution)
#Plotting in the ascending order
WaveLengthRange = 1./WaveNumberRange
WaveLengthRange = WaveLengthRange[::-1]



#Now get the assign the resolution values
Resolution = 25000

#Low resolution wavelength
Wavelength_LR, WaveNumber_LR = GetWaveNumbers(LowWavelength, HighWavelength, Resolution)

print("The range of the wavelength is given by::", Wavelength_LR)
Folder2Save = "R"+str(Resolution)

if not(os.path.exists(Folder2Save)):
    os.system("mkdir %s" %(Folder2Save))

np.savetxt(Folder2Save+"/Temperature.txt", TempRange, delimiter=",")
np.savetxt(Folder2Save+"/exp_Pressure.txt", expP_Range, delimiter=",")
np.savetxt(Folder2Save+"/WaveLength.txt", Wavelength_LR, delimiter=",")
np.savetxt(Folder2Save+"/Molecules.txt", MoleculeList, delimiter=",", fmt='%s')

BaseLocation = "DataMatrix0Sig_100cm/"
MoleculesFiles = glob.glob(BaseLocation+"*.npy")
NumMolecules = len(MoleculesFiles)
NumTempValues = len(TempRange)
NumPValues = len(expP_Range)
NumWL_Values = len(Wavelength_LR)

#Initiate a database matrix
DatabaseMatrix = np.ones((NumMolecules, NumTempValues, NumPValues, NumWL_Values), dtype=np.float32)

for MoleculeCount, Molecule in enumerate(MoleculeList):
    #Read the molecule name
    MoleculeLocation = BaseLocation+Molecule+".npy"
    print("The molecule is given by::", Molecule, ".   Now loading the data....")

    SigmaMatrix = np.load(MoleculeLocation,mmap_mode='r')
    print("Loaded the data")
    for TempCounter in range(len(TempRange)):
        for PCounter in range((len(expP_Range))):
            #Get the temperature and the pressure index...
            Sigma_HR = SigmaMatrix[TempCounter, PCounter, :][::-1]

            InterpolatedSigma  = SymplecticInterpolation(WaveLengthRange, Sigma_HR,Wavelength_LR)
            DatabaseMatrix[MoleculeCount, TempCounter, PCounter, :] = InterpolatedSigma

            plt.figure()
            plt.plot(WaveNumberRange, Sigma_HR, "ko")
            plt.plot(WaveNumber_LR, InterpolatedSigma, "r-")
            plt.show()
            print(TempCounter, PCounter)

#import the libraries
import matplotlib.pyplot as plt

#This code uses
import numpy as np
from time import time
import os
from HAPILite import CalcCrossSection


#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])    #Step size of the temperature
TempStop = float(Values[1])     #Step size of the temperature
TempStep = float(Values[2])     #Step size of the temperature

P_Start = float(Values[3])
P_Stop = float(Values[4])
P_Step = float(Values[5])                                                       #Step size of the pressure


Broadener = Values[6].replace(" ","")                                           #Broadening either self or air at this point
OmegaWidth = float(Values[7])                                                   #Consider the omegawidth -- how far the lines have to be considered
LowWavelength = float(Values[8])                                                #Shortest Wavelength coverage range
HighWavelength = float(Values[9])                                               #Longest Wavelength coverage range
WN_Resolution = float(Values[10])                                        #Resolution of the Wave Number
LineShapeProfile = Values[11].replace(" ","")                                    #Voigt profile by default
NumChunks = int(Values[12].replace(" ",""))                                      #Number of chunks for the wavenumber
MoleculeList = Values[13].split(",")                                             #Get the list of Molecular species
Cores = int(Values[14])

MoleculeList = [Item.replace(" ", "").replace("\t","") for Item in MoleculeList]


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)             #Temperature in K
P_Range = np.array([7.0,6.0,5.6,5.2,4.9,4.6,4.3,4.0,3.7,3.4,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0])-5.0

WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)            #in per cm
WaveNumberRange = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)


print("The range of temperature is given by::", TempRange)
print("The range of pressure is given by::", P_Range)



for Molecule in MoleculeList:
    print("Starting Molecule::", Molecule)
    input("Benchmark this::")
    StartTime = time()
    #initiate the saving matrix for each case
    SigmaMatrix = np.empty((len(TempRange),len(P_Range),len(WaveNumberRange)))*np.nan
    for TempCount, TempValue in enumerate(TempRange):
            for PCount, P_Value in enumerate(P_Range):
                print(TempCount, PCount)
                SigmaMatrix[TempCount, PCount, :] = CalcCrossSection(Molecule, Temp=TempValue, P = P_Value, WN_Grid=WaveNumberRange,   \
                                                    Profile="Voigt", OmegaWing=25.0, OmegaWingHW=0.0, NCORES=-1)
    #Now save the file
    np.save("DataMatrix/"+Molecule+".npy", SigmaMatrix)
    np.savetxt("DataMatrix/WN_"+str(WaveCount)+"_"+str(WaveCount+1)+".txt", WaveNumberResults)

    print("For %s, time taken is %4.2" %(Molecule, time() -StartTime))

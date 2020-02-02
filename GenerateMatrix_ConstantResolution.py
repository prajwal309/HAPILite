#import the libraries
import matplotlib.pyplot as plt

#This code uses
import numpy as np
from time import time
import os
from HAPILite import CalcCrossSection, CalcCrossSectionWithError
from lib.ReadComputeFunc import ReadData
from lib.CrossSectionFunctions import GetWaveNumbers


#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]

#Temperature and pressure is defined by the user.
TempRange = np.arange(100,901,100)
expP_Range = np.array([ 2.,1.,0.6,0.2,-0.1,-0.4,-0.7,-1.,-1.3, \
                         -1.6,-2.,-3.,-4.,-5.,-6.,-7.,-8.])

OmegaWidth = 100.0                                                              #Consider the omegawidth -- how far the lines have to be considered
LowWavelength = float(Values[8])                                                #Shortest Wavelength coverage range
HighWavelength = float(Values[9])                                               #Longest Wavelength coverage range
LineShapeProfile = Values[11].replace(" ","")                                   #Voigt profile by default
MoleculeList = Values[12].split(",")                                            #Get the list of Molecular species
Cores = int(Values[13])
Error = "1SIG"#Values[14].replace("\t","")

SaveFolder = "DataMatrix"+Error.replace("-","Neg").replace(" ","")

if not(os.path.exists(SaveFolder)):
    os.system("mkdir %s" %SaveFolder)


MoleculeList = [Item.replace(" ", "").replace("\t","") for Item in MoleculeList]
MoleculeList = np.array([Item.replace(" ", "").replace("\t","") for Item in MoleculeList])


#Now get the assign the resolution values
Resolution = 100000

#Low resolution wavelength
WavelengthRange, WaveNumberRange = GetWaveNumbers(LowWavelength, HighWavelength, Resolution)

for Molecule in MoleculeList:
    print("\n\n Starting Molecule::", Molecule)
    StartTime = time()
    Database = ReadData(Molecule, Location="data/")

    #initiate the saving matrix for each case
    SigmaMatrix = np.zeros((len(TempRange),len(expP_Range),len(WaveNumberRange)),dtype=np.float32)
    print("The shape of the Sigma Matrix is given by:", np.shape(SigmaMatrix))

    for TempCount, TempValue in enumerate(TempRange):
            for PCount, expPValue in enumerate(expP_Range):

                P_Value = 10**expPValue
                print("-"*15)
                print("Temperature:", TempValue)
                print("Pressure:", P_Value)

                if "0SIG" in Error.upper():
                    SigmaMatrix[TempCount, PCount, :] = CalcCrossSection(Database, Temp=TempValue, P = P_Value, WN_Grid=WaveNumberRange,   \
                                                        Profile=LineShapeProfile, OmegaWing=100.0, OmegaWingHW=0.0, NCORES=Cores)[::-1]

                elif "1SIG" in Error.upper() or "2SIG" in Error.upper():
                    SigmaMatrix[TempCount, PCount, :] = CalcCrossSectionWithError(Database, Temp=TempValue, P = P_Value, WN_Grid=WaveNumberRange,   \
                                                        Profile=LineShapeProfile, OmegaWing=100.0, OmegaWingHW=0.0, NCORES=Cores, Err=Error)[::-1]

                else:
                    raise Exception("Error in GenerateMatrix. The error has to be 0SIG, 1SIG/-1SIG or 2SIG/-2SIG ")

    #Now save the file


    np.save(SaveFolder+"/"+Molecule+".npy", SigmaMatrix)
    print("For %s, time taken is %5.2f" %(Molecule, time() -StartTime))

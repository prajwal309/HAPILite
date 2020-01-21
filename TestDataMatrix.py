import numpy as np
import time
import hapi
from HAPILite import CalcCrossSection
from lib.ReadComputeFunc import ReadData
import matplotlib.pyplot as plt
import glob


Molecule = "CO"
TempValue = 800.0
P_Value = 10.0
Env= {'T': TempValue, 'p': P_Value}


#
#
#Read the data from the location

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

print("The omega value::",OmegaWidth)
input("Wait here...")

WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)               #in per cm
WaveNumberRange = np.arange(WaveNumberStart, WaveNumberStop, WN_Resolution)


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)                          #Temperature in K
P_Range = 10.0**np.arange(expP_Start, expP_Stop-expP_Step, -expP_Step)                #Pressure in log(P) atm

#Calculate from the HAPI
hapi.db_begin("TestData")
nu_hapi_Voigt, abs_hapi_Voigt = hapi.absorptionCoefficient_Voigt(SourceTables=Molecule,OmegaGrid=WaveNumberRange,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', OmegaWing=OmegaWidth, OmegaWingHW=0.0, LineShift=False)#, OmegaStep=0.01)


#Get the Index
TemperatureIndex = np.argmin(np.abs(TempRange - TempValue))
PressureIndex = np.argmin(np.abs(P_Range - P_Value))

#Read the data
FileLocation = glob.glob("DataMatrix0Sig_100cm/%s*" %(Molecule))[0]
Data = np.load(FileLocation, mmap_mode='r')

CorrespondingData = Data[TemperatureIndex, PressureIndex, :]

fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)
ax0.plot(nu_hapi_Voigt, abs_hapi_Voigt, "k.", label="Old HAPI" )
ax0.plot(WaveNumberRange, CorrespondingData, "r-",linewidth=1.0, label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")
ax0.set_ylabel("Cross-Section (per cm)")
ax0.legend(loc=1)
#ax1.plot(nu_hapi_Voigt, abs_hapi_Voigt - CorrespondingData, "r-" )
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("Voigt ---   %s (Temp: %s Pressure: %s)" %(Molecule, str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
plt.savefig("CrossSectionTest_%s_Voigt.png" %(Molecule))
plt.close()

import numpy as np
import matplotlib.pyplot as plt
from lib.ReadComputeFunc import ReadData



#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])    #Step size of the temperature
TempStop = float(Values[1])     #Step size of the temperature
TempStep = 25                   #Step size of the temperature

expP_Start = float(Values[3])   #The largest log10(pressure) in atm
expP_Stop = float(Values[4])    #The smallest log10(pressure) in atm
expP_Step = 0.25    #Step size of the pressure

LowWavelength = float(Values[8])                                                #Shortest Wavelength coverage range
HighWavelength = float(Values[9])                                               #Longest Wavelength coverage range
WN_Resolution = float(Values[10])

WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)               #in per cm
WaveNumberRange = np.arange(WaveNumberStart, WaveNumberStop, WN_Resolution)


print("The wavenumber start is::", WaveNumberStart)
print("The wavenumber start is::", WaveNumberStop)
print("The length of the wavenumber range is::", len(WaveNumberRange))


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)                    #Temperature in K
expP_Range = np.arange(expP_Start, expP_Stop-expP_Step, -expP_Step)             #Pressure in log(P) atm


Database = ReadData("H2O")




#Read the cross -section
MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaSelf, TempRatioPower, ErrorArray = Database


#Total sum of the wavenumber is given by
SelectIndex = np.logical_and(LineCenterDB>WaveNumberStart, LineCenterDB<WaveNumberStop)

SelectedCenterDB = LineCenterDB[SelectIndex]
SelectedIntensityDB = LineIntensityDB[SelectIndex]

TotalIntensity = np.sum(SelectedIntensityDB)


ErrorMatrix1 = np.zeros((len(TempRange), len(expP_Range)))
ErrorMatrix2 = np.zeros((len(TempRange), len(expP_Range)))

print("The value of error matrix is given by")

Database1 = np.load("DataMatrix0Sig_25cm/H2O.npy", mmap_mode='r' )
Database2 = np.load("DataMatrix0Sig_100cm/H2O.npy", mmap_mode='r')

cBolts = 1.380648813E-16 # erg/K, CGS
cc = 2.99792458e10 # cm/s, CGS
hh = 6.626196e-27 # erg*s, CGS

for TCounter, Temp in enumerate(TempRange):
    for PCounter, expP in enumerate(expP_Range):

        PValue = 10.**expP
        if (PValue-1.0)<0.1:
            print("The value of pressure is::")
            factor = (PValue/9.869233e-7)/(cBolts*Temp) #
            print("The value of factor is::", factor)


            AbsCrossSection1 = Database1[TCounter, PCounter,:]/factor
            AbsCrossSection2 = Database2[TCounter, PCounter,:]/factor

            plt.figure()
            plt.plot(WaveNumberRange, AbsCrossSection1, "k.")
            plt.plot(SelectedCenterDB, SelectedIntensityDB, "g.")
            plt.show()

            #Integration value of the database
            Integration1 = np.trapz(AbsCrossSection1, WaveNumberRange)*296./Temp
            Integration2 = np.trapz(AbsCrossSection2, WaveNumberRange)*296./Temp

            RelativeError1 = (TotalIntensity-Integration1)/TotalIntensity*100.0
            RelativeError2 = (TotalIntensity-Integration2)/TotalIntensity*100.0

            print("*"*25)
            print("The value of relative error 25 cm is::", round(RelativeError1,5))
            print("The value of relative error 100 cm is::", round(RelativeError2,5))
            print("*"*25)

            ErrorMatrix1[TCounter, PCounter] = RelativeError1
            ErrorMatrix2[TCounter, PCounter] = RelativeError2



plt.figure()
plt.imshow(ErrorMatrix1, origin="lower")
plt.show()

import glob
import numpy as np
from CrossSectionFunctions import GetWaveNumbers, SymplecticInterpolation
import matplotlib.pyplot as plt
import os
#from numba import jit, njit


#parse the parameters.dat
Data = [f.split(":") for f in open("Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])                                                     #Step size of the temperature
TempStop = float(Values[1])                                                     #Step size of the temperature
TempStep = float(Values[2])                                                     #Step size of the temperature

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

MoleculeList = [Item.replace(" ", "").replace("\t","") for Item in MoleculeList]


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)             #Temperature in K
#P_Range = np.arange(P_Start, P_Stop+P_Step, P_Step)        #log10(Pressure) in atm
P_Range = np.array([7.0,6.0,5.6,5.2,4.9,4.6,4.3,4.0,3.7,3.4,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0])-5.0




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
np.savetxt(Folder2Save+"/Pressure.txt", P_Range)
np.savetxt(Folder2Save+"/WaveLength.txt", Wavelength_LR)


#converting to nm
Wavelength_LR*=1e7


LengthWaveNumber = len(np.arange(WaveNumberStart, WaveNumberStop+WN_Resolution, WN_Resolution))

WaveNumberRanges = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)


#os.system("rm HereBeData/*")
#mkdir("1000")
print("Removing the data...")

for Molecule in MoleculeList:
    #There is only

    FileList = np.array(FileList)[np.argsort(Index)]
    StartPoint = WaveNumberStart

    print("The files after arranged file list is::")
    print(FileList)

    #Generate a file per molecule
    X,Y,Z = len(TempRange), len(P_Range), len(Wavelength_LR)
    print("The value of X is given::", X)
    print("The value of Y is given::", Y)

    input("Crash here")

    MoleculeMatrix = np.zeros((X,Y,Z))

    for FileCount,File in enumerate(FileList):
        WaveNumberIndex = np.logical_and(WaveNumber_LR>WaveNumberRanges[FileCount], WaveNumber_LR<WaveNumberRanges[FileCount+1])
        SelectedWaveNumber = WaveNumber_LR[WaveNumberIndex]

        Data = np.load(File)

        #Now generate the right cross section
        m, n, k = np.shape(Data)

        print("The name of the filecount is::", FileCount)

        for OuterCounter in range(m):
            for InnerCounter in range(n):
                GenCrossSection = Data[OuterCounter,InnerCounter,:]

                TemperatureValue = TempRange[OuterCounter]
                PressureValue = round(P_Range[InnerCounter],3)

                Selected_WN = np.linspace(WaveNumberRanges[FileCount], WaveNumberRanges[FileCount+1], len(GenCrossSection))

                InterpolatedValues = SymplecticInterpolation(Selected_WN,GenCrossSection, SelectedWaveNumber)

                #Convert WaveNumbers to Wavelengths
                Selected_HR_Wavelength = 1./Selected_WN*1e7
                Selected_LR_Wavelength = 1./SelectedWaveNumber*1e7

                MinWavelength  = min(Selected_HR_Wavelength)
                MaxWavelength = max(Selected_HR_Wavelength)

                #Reverse the sigma when assigning the value
                AssignIndexStart = np.argmin(np.abs(MinWavelength-Wavelength_LR))
                AssignIndexStop = AssignIndexStart + len(InterpolatedValues)

                #Check if the index match
                #if np.argmin(np.abs(MaxWavelength-Wavelength_LR))

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




    SaveLocation = Folder2Save+"/"+Molecule+".npy"
    print(SaveLocation)

    #Now save the data
    #Findif any nan
    NanIndex = np.sum(np.isnan(MoleculeMatrix))
    print("The nanindex, for ", Molecule, " is given by ", NanIndex)
    np.save(SaveLocation, MoleculeMatrix)
    #input("Wait here...")

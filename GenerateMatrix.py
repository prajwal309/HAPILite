#import the libraries
import matplotlib.pyplot as plt

#This code uses
import numpy as np
from itertools import product
from time import time
import numba as nb
import multiprocessing as mp
import os



#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
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
Cores = int(Values[14])

MoleculeList = [Item.replace(" ", "").replace("\t","") for Item in MoleculeList]


TempRange = np.arange(TempStart,TempStop+TempStep, TempStep)             #Temperature in K
P_Range = np.array([7.0,6.0,5.6,5.2,4.9,4.6,4.3,4.0,3.7,3.4,3.0,2.0,1.0,0.0,-1.0,-2.0,-3.0])-5.0

WaveNumberStart = 1./(HighWavelength*1.e-7)            #in per cm
WaveNumberStop= 1./(LowWavelength*1.e-7)            #in per cm
WaveNumberRange = np.linspace(WaveNumberStart, WaveNumberStop, NumChunks+1)

for WaveCount in range(NumChunks):
    WN_Start = round(WaveNumberRange[WaveCount],4)
    WN_Stop = round(WaveNumberRange[WaveCount+1],4)
    if WaveCount == 0:
        #Single Length to rule them all
        Length = len(np.arange(WN_Start, WN_Stop, WN_Resolution))

    SigmaMatrix = np.empty((len(TempRange),len(P_Range),Length))*np.nan

    StartTime = time()

    for Molecule in MoleculeList:

        #copy data to a temperorary folder
        if not(os.path.exists("./TempData")):
            os.system("mkdir ./TempData")
        os.system("rm ./TempData/*")
        os.system("cp ./data/%s.* ./TempData" %(Molecule))
        TABLE_NAME = os.getcwd()+"/TempData/"+Molecule
        hapi.storage2cache(TABLE_NAME)
        hapi.VARIABLES['CPF'] = hapi.hum1_wei


        if "self" in Broadener:
            DILUENT = [('self',1.0)]
        elif "air" in Broadener:
            DILUENT = [('air',1.0)]
        else:
            print("Other broadeners have not been implemented yet.")
            assert 1==2


        MOLEC_ID,LOCAL_ISO_ID = hapi.getColumns(TABLE_NAME,['molec_id','local_iso_id'])
        NLINES = len(MOLEC_ID)
        ISOS = GET_ISOS_DEFAULT_ABUN(NLINES,MOLEC_ID,LOCAL_ISO_ID)
        MOLEC_ID,LOCAL_ISO_ID = hapi.getColumns(TABLE_NAME,['molec_id','local_iso_id'])
        NCORES = 1
        OMEGA_STEP = WN_Resolution

        OMEGA_RANGE = np.array([WN_Start, WN_Stop])
        OMEGAWING = 0.0
        OMEGAWINGHW = OmegaWidth
        # 1 - Voigt, 2 - Lorentz, 3 - Doppler
        if "voigt" in LineShapeProfile.lower():
            PROFILE_NUM = 1
        elif "lorentz" in LineShapeProfile.lower():
            PROFILE_NUM = 2
        elif "doppler" in LineShapeProfile.lower():
            PROFILE_NUM = 3
        else:
            print("Only profile number 1,2 and 3 are available.")
            assert 1==2


        NCORES_ABS = 1
        if Cores<1:
            NUM_CORES  = mp.cpu_count()
        else:
            NUM_CORES = Cores
        #Database for storing the value
        Product_TP = product(range(len(TempRange)), range(len(P_Range)))
        for counter in range(int(len(TempRange)*len(P_Range)/NUM_CORES)+1):
            print("The value of counter is:", counter)
            CPU_Pool = mp.Pool(NUM_CORES)
            #List to temporily store values
            KeyWordsList = []           #Keyword arguments
            TempCountList = []          #Temperature Count List
            PCountList = []             #Pressure Count List

            #Start as many process as number of cores
            #Just generate the relevant keywords
            for MPICount in range(NUM_CORES):
                try:
                    TempCount, PCount = next(Product_TP)
                except:
                    pass

                TempCountList.append(TempCount)
                PCountList.append(PCount)

                TemperatureValue = TempRange[TempCount]
                PressureValue = 10**(P_Range[PCount])               #Converting to atm

                #TempKEYWORDS = {"OmegaWingHW":OMEGAWINGHW,
                #"OmegaRange":OMEGA_RANGE,"OmegaStep":OMEGA_STEP,
                #"reflect":False,"T":TemperatureValue,"p":PressureValue,"profile":PROFILE_NUM,"NCORES":1}

                TempKEYWORDS = {"OmegaWing":OMEGAWING,"OmegaWingHW":OMEGAWINGHW,
                "OmegaRange":OMEGA_RANGE,"OmegaStep":OMEGA_STEP,
                "reflect":False,"T":TemperatureValue,"p":PressureValue,"profile":PROFILE_NUM,"NCORES":1}
                KeyWordsList.append(TempKEYWORDS)

            #Append the tasks
            Tasks = [CPU_Pool.apply_async(ABSCOEF_FAST, (NLINES,TABLE_NAME,ISOS,DILUENT), KeyWordsList[x]) for x in range(MPICount+1)]

            #Following function wraps up the CPU functions and completes the processes
            CPU_Pool.close()
            CPU_Pool.join()

            for AssignCount, task in enumerate(Tasks):


                #print("The AssignCount is::", AssignCount)


                Results = np.array(task.get())
                Sigma = Results[1,:]
                if len(Sigma)>Length:
                    print("Error is the length::",len(Sigma))
                    Sigma = Sigma[:Length]
                SigmaMatrix[TempCountList[AssignCount],PCountList[AssignCount],:] = Sigma
                WaveNumberResults = Results[0,:]


                TValue = TempRange[TempCountList[AssignCount]]
                PValue = round(10.0**P_Range[PCountList[AssignCount]],2)



        if not(os.path.exists("DataMatrix")):
            os.system("mkdir DataMatrix")

        #Now save the file
        np.save("DataMatrix/"+Molecule+"_"+str(WaveCount)+"_"+str(WaveCount+1)+".npy", SigmaMatrix)
        np.savetxt("DataMatrix/WN_"+str(WaveCount)+"_"+str(WaveCount+1)+".txt", WaveNumberResults)
    StopTime = time()
    print("The time taken for ", Molecule, "is ", StopTime - StartTime)

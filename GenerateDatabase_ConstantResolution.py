import glob
import numpy as np
from lib.CrossSectionFunctions import GetWaveNumbers
import matplotlib.pyplot as plt
import os
import itertools
import time
import multiprocessing as mp

#Temperature and pressure is defined by the user.
TempRange = np.arange(100,901,100)
expP_Range = np.array([ 2.,1.,0.6,0.2,-0.1,-0.4,-0.7,-1.,-1.3, \
                         -1.6,-2.,-3.,-4.,-5.,-6.,-7.,-8.])

#Now get the assign the resolution values
Resolution = 100000

#Low resolution wavelength
WavelengthRange, WaveNumberRange = GetWaveNumbers(300, 30000, Resolution)

BaseLocation = "DataMatrix1SIG"
Folder2Save = "New_R100000_1SIG"

if not(os.path.exists(Folder2Save)):
    os.system("mkdir %s" %Folder2Save)

FileList = glob.glob(BaseLocation+"/*.npy")
print("The list of the files is given by::", FileList)
input("Wait here...")
MoleculeList = [Item.split("/")[-1][:-4] for Item in FileList]


np.savetxt(Folder2Save+"/Temperature.txt", TempRange, delimiter=",")
np.savetxt(Folder2Save+"/exp_Pressure.txt", expP_Range, delimiter=",")
np.savetxt(Folder2Save+"/WaveLength.txt", WavelengthRange, delimiter=",")
np.savetxt(Folder2Save+"/Molecules.txt", MoleculeList, delimiter=",", fmt='%s')

CSMatrix = np.zeros((len(TempRange), len(expP_Range), len(MoleculeList), len(WavelengthRange)),dtype=np.float32)

print("The shape of the CSMatrix is::", np.shape(CSMatrix))

for MolPos, File in enumerate(FileList):
    CurrentData = np.load(File, mmap_mode='r')
    CSMatrix[:,:,MolPos,:] = CurrentData

#input("becuase Forgot to invert the row")
#for MolPos, File in enumerate(FileList):
#    CurrentData = np.load(File, mmap_mode='r')
#    for TempCount in range(len(TempRange)):
#        for PCount in range(len(expP_Range)):
#            RowData = CurrentData[TempCount, PCount,:][::-1]
#            CSMatrix[TempCount, PCount,MolPos,:] = RowData

np.save(Folder2Save+"/DataBase_%d.npy" %Resolution, CSMatrix)

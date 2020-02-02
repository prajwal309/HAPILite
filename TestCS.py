import numpy as np
import matplotlib.pyplot as plt

#Read data
#Location1 = "New_R100000_1SIG"
Location1 = "/home/prajwal/Desktop/HITRAN/SpectroscopyHITRAN/PythonAdaptation/DatabaseHL/New_R100000"
Location2 = "/home/prajwal/Desktop/HITRAN/SpectroscopyHITRAN/PythonAdaptation/DatabaseHL/R10000"

Wavelength1 = np.loadtxt(Location1+"/"+"WaveLength.txt")
Wavelength2 = np.loadtxt(Location2+"/"+"WaveLength.txt")

TempVect1 = np.loadtxt(Location1+"/"+"Temperature.txt")
TempVect2 = np.loadtxt(Location2+"/"+"Temperature.txt")

Pressure1 = np.loadtxt(Location1+"/"+"exp_Pressure.txt")
Pressure2 = np.loadtxt(Location2+"/"+"exp_Pressure.txt")

Molecule1 = np.loadtxt(Location1+"/"+"Molecules.txt", dtype=np.str)
Molecule2 = np.loadtxt(Location2+"/"+"Molecules.txt", dtype=np.str)

TValue = 800
exp_PValue = 2.0
MoleculeValue = "N2"

TIndex1 = np.where(TempVect1==TValue)[0][0]
PIndex1 = np.where(Pressure1==exp_PValue)[0][0]
MoleculeIndex1 = np.where(Molecule1 == MoleculeValue)[0][0]

TIndex2 = np.where(TempVect2==TValue)[0][0]
PIndex2 = np.where(Pressure2==exp_PValue)[0][0]
MoleculeIndex2 = np.where(Molecule2 == MoleculeValue)[0][0]

print("Loading data...")
Database1 = np.load(Location1+"/Database.npy", mmap_mode='r')
Database2 = np.load(Location2+"/Database.npy", mmap_mode='r')

Sigma1 = Database1[TIndex1, PIndex1, MoleculeIndex1, :]
Sigma2 = Database2[TIndex2, PIndex2, MoleculeIndex2, :]




plt.figure(figsize=(12,8))
plt.plot(Wavelength1*1e7, Sigma1/np.max(Sigma1), "ko", label="New Cross-Section")
plt.plot(Wavelength2*1e7, Sigma2/np.max(Sigma2), "r-", label="Old Cross-Section")

TitleText = Molecule1[MoleculeIndex1]+":"+"  T:"+str(TempVect1[TIndex1])+"  P:"+str(Pressure1[PIndex1])

plt.title(TitleText)
plt.xscale('log')
plt.xlabel("Wavelength (nm)")

plt.legend(loc=1)
plt.savefig("ComparisonFigure.png")
plt.show()

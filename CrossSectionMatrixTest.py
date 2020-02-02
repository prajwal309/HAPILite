import numpy as np
import time
import matplotlib.pyplot as plt
from HAPILite import CalcCrossSection, CalcCrossSectionWithError
from lib.ReadComputeFunc import ReadData
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter

import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})




Molecule = "H2O"
TempValue = 600.0
P_Value = 1.0

OmegaWingValue = 100.0
OmegaRangeValue = [1./30000.*1e7, 1./300.*1e7]
print("The omega range values is given by::", OmegaRangeValue)


WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.01, 0.005)


StartTime = time.time()

Database = ReadData(Molecule, Location="data/")
MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, \
LowerStateEnergyDB, GammaSelf, TempRatioPower, ErrorArray = Database

SelectIndex = np.logical_and(LineCenterDB>OmegaRangeValue[0], LineCenterDB<OmegaRangeValue[1])

LineCenterDB = LineCenterDB[SelectIndex]
LineIntensityDB = LineIntensityDB[SelectIndex]

ch = np.exp(-const_R*LowerStateEnergy/Temp)*(1-np.exp(-const_R*LineCenterDB/Temp))
zn = np.exp(-const_R*LowerStateEnergy/Tref)*(1-np.exp(-const_R*LineCenterDB/Tref))
LineIntensity = LineIntensityDB[i]*SigmaTref/SigmaT*ch/zn

TotalIntensity = np.sum(LineIntensityDB)


CrossSection =  CalcCrossSectionWithIntensity(Database,Temp=TempValue, P = 1.0,\
                         WN_Grid=WaveNumber, Profile="Doppler", OmegaWing=OmegaWingValue,\
                         OmegaWingHW=0.0, NCORES=-1, Err="0SIG")


Area = np.trapz(CrossSection, WaveNumber)

#Generate the area

print("The total value of intensity is::", np.sum(TotalIntensity))
print("The value of the area is given by::", Area)
input("Wait here...")


plt.figure()
plt.plot(WaveNumber, CrossSection, "k.")
plt.savefig("DeleteMe.png")
plt.close('all')

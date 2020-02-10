
import numpy as np
import matplotlib.pyplot as plt
from HAPILite import CalcCrossSection
from lib.ReadComputeFunc import ReadData
from scipy.stats import pearsonr


import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})

Molecule = "CH4"
Database = ReadData(Molecule, Location="data/")

MoleculeNumberDB, IsoNumberDB, LineCenterDB, LineIntensityDB, LowerStateEnergyDB, GammaAir, GammaSelf, TempRatioPower, ErrorArray = Database


CorrValue1 = pearsonr(GammaAir, GammaSelf)
CorrValue2 = pearsonr(GammaAir, GammaSelf)

plt.figure(figsize=(15,8))
plt.subplot(211)
plt.ylabel("Pressure Broadening Coefficient")
plt.plot(LineIntensityDB, GammaSelf, "ko")
plt.xlabel("Line Intensity")
plt.ylabel("Pressure Broadening")

plt.xscale('log')
plt.subplot(212)
plt.plot(GammaAir, GammaSelf, "ko")
plt.suptitle(Molecule)
plt.xlabel("Air Broadening")
plt.ylabel("Self Broadening")
plt.tight_layout()
plt.savefig(Molecule+"BroadeningCoefficient.png")

plt.close()

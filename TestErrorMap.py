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
OmegaRangeValue = [2000, 2060]

WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.01, 0.01)


StartTime = time.time()

Database = ReadData(Molecule, Location="data/")


CrossSectionVoigt0Sig =  CalcCrossSectionWithError(Database,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=-1, Err="0SIG")

CrossSectionVoigt1SigPos =  CalcCrossSectionWithError(Database,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=-1, Err="1SIG")

CrossSectionVoigt1SigNeg =  CalcCrossSectionWithError(Database,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=-1, Err="-1SIG")


fig, ax = plt.subplots(figsize=(6,6), dpi=200)
ax.plot(WaveNumber, CrossSectionVoigt0Sig, "k-", lw=1.5, label = "Standard Cross-Section")
ax.fill_between(WaveNumber, CrossSectionVoigt1SigPos, CrossSectionVoigt1SigNeg, color="red", alpha=0.90, label="$\\pm$1$\\sigma$")
ax.set_ylabel("Cross-Section [cm$^{-1}$/mol cm$^{-2}$]", fontsize=20)
ax.set_ylim([1.20e-21, 1.20e-19])
ax.yaxis.set_major_formatter(FormatStrFormatter('%3.2e'))
ax.set_xlabel("Wavenumber (cm$^{-1}$)", fontsize=20)
#ax.set_xlim([1660.0, 1750.0])
#ax.set_xlim([1600.0, 1690.0])
plt.tick_params(which="both", direction="in")
#plt.title("Voigt   %s (Temp: %s Pressure: %s)" %(Molecule, str(int(TempValue)),str(round(P_Value,3))))
plt.legend()
plt.tight_layout()
plt.savefig("ErrorFigure.png")
plt.savefig("ErrorFigure.pdf", format="pdf")
plt.show()
plt.close()

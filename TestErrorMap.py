import numpy as np
import time
import matplotlib.pyplot as plt
from HAPILite import CalcCrossSection, CalcCrossSectionWithError
from lib.ReadComputeFunc import ReadData

#Molecule = "N2"
#TempValue = 900.0
#P_Value = 10.0

#OmegaWingValue = 25.0
#OmegaRangeValue = [2301.6, 2302.6]



Molecule = "H2O"
TempValue = 900.0
P_Value = 10.0

OmegaWingValue = 100.0
#OmegaRangeValue = [2049.5, 2070.1]
OmegaRangeValue = [0, 10000]

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


plt.figure(figsize=(12,6))
plt.plot(WaveNumber, CrossSectionVoigt0Sig, "k-", alpha=0.75)
plt.fill_between(WaveNumber, CrossSectionVoigt1SigPos, CrossSectionVoigt1SigNeg, color="red", alpha=0.90)
plt.tick_params(which="both", direction="in")
plt.ylabel("Cross-Section")
plt.xlabel("Wavenumber (cm$^{-1}$)")
plt.title("Voigt   %s (Temp: %s Pressure: %s)" %(Molecule, str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
plt.savefig("ErrorFigure.png")
plt.show()
plt.close('all')

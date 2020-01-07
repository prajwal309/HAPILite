import numpy as np
import time
\
from HAPILite import CalcCrossSection, CalcCrossSectionWithError
import matplotlib.pyplot as plt


Molecule = "N2"
TempValue = 900.0
P_Value = 1.0

OmegaWingValue = 25.0
OmegaRangeValue = [333.333,33333.33333]
WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.01, 0.01)





print("\n"*10)
print("*"*25)
StartTime = time.time()

CrossSectionVoigt0Sig =  CalcCrossSectionWithError(Molecule,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=1, Err="0SIG")

CrossSectionVoigt1SigPos =  CalcCrossSectionWithError(Molecule,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=1, Err="1SIG")

CrossSectionVoigt1SigNeg =  CalcCrossSectionWithError(Molecule,Temp=TempValue, P = P_Value,\
                     WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue,\
                     OmegaWingHW=0.0, NCORES=1, Err="-1SIG")


plt.figure()
plt.plot(WaveNumber, CrossSectionVoigt0Sig, "k.")
plt.plot(WaveNumber, CrossSectionVoigt1SigPos, "r-")
plt.plot(WaveNumber, CrossSectionVoigt1SigNeg, "r-")
plt.show()

#author:Prajwal Niraula
#insitution: MIT


import matplotlib.pyplot as plt
import time
import numpy as np

try:
    from HAPILite import CalcCrossSection
except:
    from ..HAPILite import CalcCrossSection

WaveNumber = np.arange(0,10000,0.001)

StartTime = time.time()
CrossSection =  CalcCrossSection("CO2",Temp=1000.0,WN_Grid=WaveNumber, Profile="Doppler", NCORES=-1)
StopTime = time.time()

print("The time take to calculate the cross-section is %4.3f" %(StopTime - StartTime))


plt.figure()
plt.plot(WaveNumber, CrossSection, "k-")
plt.title("L-HAPI")
plt.show()

import numpy as np
import time
import hapi
from HAPILite import CalcCrossSection
import matplotlib.pyplot as plt


TempValue = 300.0
P_Value = 1.0
OmegaRangeValue = [2196.93, 2197.02]
WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1], 0.0005)
Env= {'T': TempValue, 'p': P_Value}


hapi.db_begin('data')
StartTime = time.time()
nu_hapi, abs_hapi = hapi.absorptionCoefficient_Voigt(SourceTables='CO2',OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', WavenumberWingHW=75.0)#, OmegaStep=0.01)
TimeTakenHapi = time.time() - StartTime
print("The time taken for HAPI is::",TimeTakenHapi)

StartTime = time.time()
CrossSection =  CalcCrossSection("CO2",Temp=1000.0,WN_Grid=WaveNumber, Profile="Voigt", NCORES=-1)
TimeTakenHAPILite = time.time() - StartTime
print("The time taken by HAPI Lite is::", TimeTakenHAPILite)

print("The time difference is::", TimeTakenHapi - TimeTakenHAPILite)

print("The average multiplicative factor is::", np.mean(abs_hapi/CrossSection))

fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)

ax0.plot(nu_hapi, abs_hapi, "ko-", label="Old HAPI" )
ax0.plot(WaveNumber, CrossSection, "ro:", label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")
ax0.set_ylabel("Cross-Section (per cm)")

ax1.plot(nu_hapi, abs_hapi/CrossSection, "r-" )
ax1.set_yscale("log")
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("CO2 (Temp: %s Pressure: %s)" %(str(int(TempValue)),str(round(P_Value,3))))
#plt.tight_layout()
plt.savefig("BenchmarkTestCO2.png")
plt.show()

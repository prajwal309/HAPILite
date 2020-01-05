import numpy as np
import time
import hapi
from HAPILite import CalcCrossSection
import matplotlib.pyplot as plt


TempValue = 1000.0
P_Value = 1.0
OmegaRangeValue = [0, 20000]
WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1], 0.01)
Env= {'T': TempValue, 'p': P_Value}


hapi.db_begin('data')
nu_hapi_Doppler, abs_hapi_Doppler = hapi.absorptionCoefficient_Doppler(SourceTables='CO2',OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', WavenumberWingHW=75.0, LineShift=False)#, OmegaStep=0.01)
StartTime = time.time()
nu_hapi_Voigt, abs_hapi_Voigt = hapi.absorptionCoefficient_Voigt(SourceTables='CO2',OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', WavenumberWingHW=75.0, LineShift=False)#, OmegaStep=0.01)
TimeTakenHapi = time.time() - StartTime
print("The time taken for HAPI is::",TimeTakenHapi)

CrossSectionDoppler =  CalcCrossSection("CO2",Temp=1000.0,WN_Grid=WaveNumber, Profile="Doppler", OmegaWingHW=75.0, NCORES=-1)
StartTime = time.time()
CrossSectionVoigt =  CalcCrossSection("CO2",Temp=1000.0,WN_Grid=WaveNumber, Profile="Voigt", OmegaWingHW=75.0, NCORES=-1)
TimeTakenHAPILite = time.time() - StartTime
print("The time taken by HAPI Lite is::", TimeTakenHAPILite)



print("The time difference is::", TimeTakenHapi - TimeTakenHAPILite)

fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)

ax0.plot(nu_hapi_Doppler, abs_hapi_Doppler, "ko-", label="Old HAPI" )
ax0.plot(WaveNumber, CrossSectionDoppler, "ro:", label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")
ax0.set_ylabel("Cross-Section (per cm)")
ax0.legend(loc=1)
ax1.plot(nu_hapi_Doppler, abs_hapi_Doppler - CrossSectionDoppler, "r-" )
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("Doppler   CO2 (Temp: %s Pressure: %s)" %(str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
plt.savefig("BenchmarkTestCO_Doppler.png")
plt.close()


fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)

ax0.plot(nu_hapi_Voigt, abs_hapi_Voigt, "ko-", label="Old HAPI" )
ax0.plot(WaveNumber, CrossSectionVoigt, "ro:", label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")

ax0.set_ylabel("Cross-Section (per cm)")
ax0.legend(loc=1)
ax1.plot(nu_hapi_Voigt, abs_hapi_Voigt - CrossSectionVoigt, "r-" )
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("Voigt ---   CO2 (Temp: %s Pressure: %s)" %(str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
plt.savefig("BenchmarkTestCO_Voigt.png")
plt.close()

import numpy as np
import time
import hapi
from HAPILite import CalcCrossSection
import matplotlib.pyplot as plt


Molecule = "H2O"
TempValue = 500.0
P_Value = 0.1

#TempValue = 1000.0
#P_Value = 10.0
OmegaWingValue = 0.1
OmegaRangeValue = [333.333,33333.33333]
WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.001, 0.001)
Env= {'T': TempValue, 'p': P_Value}


hapi.db_begin('TestData')
nu_hapi_Doppler, abs_hapi_Doppler = hapi.absorptionCoefficient_Doppler(SourceTables=Molecule,OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', OmegaWing=OmegaWingValue, OmegaWingHW=OmegaWingValue, LineShift=False)#, OmegaStep=0.01)
StartTime = time.time()
nu_hapi_Voigt, abs_hapi_Voigt = hapi.absorptionCoefficient_Voigt(SourceTables=Molecule,OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', OmegaWing=OmegaWingValue, OmegaWingHW=OmegaWingValue, LineShift=False)#, OmegaStep=0.01)
TimeTakenHapi = time.time() - StartTime
print("The time taken for HAPI is::",TimeTakenHapi)

CrossSectionDoppler =  CalcCrossSection(Molecule,Temp=TempValue,P = P_Value, WN_Grid=WaveNumber, Profile="Doppler", OmegaWing=OmegaWingValue, OmegaWingHW=0.0, NCORES=-1)


print("\n"*10)
print("*"*25)
StartTime = time.time()
CrossSectionVoigt =  CalcCrossSection(Molecule,Temp=TempValue, P = P_Value, WN_Grid=WaveNumber, Profile="Voigt", OmegaWing=OmegaWingValue, OmegaWingHW=0.0, NCORES=-1)
TimeTakenHAPILite = time.time() - StartTime
print("The time taken by HAPI Lite is::", TimeTakenHAPILite)



print("The time difference is::", TimeTakenHapi - TimeTakenHAPILite)

fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)
ax0.plot(nu_hapi_Doppler, abs_hapi_Doppler, "k.", label="Old HAPI" )
ax0.plot(WaveNumber, CrossSectionDoppler, "r-", linewidth=1.0, label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")
ax0.set_ylabel("Cross-Section (per cm)")
ax0.legend(loc=1)
ax1.plot(nu_hapi_Doppler, abs_hapi_Doppler - CrossSectionDoppler, "r-" )
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("Doppler   %s (Temp: %s Pressure: %s)" %(Molecule, str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
plt.savefig("BenchmarkTest_%s_Doppler.png" %(Molecule))
plt.close()


fig, (ax0, ax1) = plt.subplots(figsize=(14,6), ncols=1, nrows=2, sharex=True)

ax0.plot(nu_hapi_Voigt, abs_hapi_Voigt, "k.", label="Old HAPI" )
ax0.plot(WaveNumber, CrossSectionVoigt, "r-",linewidth=1.0, label="HAPI Lite" )
ax0.set_xlabel("WaveNumber (per cm)")

ax0.set_ylabel("Cross-Section (per cm)")
ax0.legend(loc=1)
ax1.plot(nu_hapi_Voigt, abs_hapi_Voigt - CrossSectionVoigt, "r-" )
ax1.set_xlabel("WaveNumber (per cm)")
ax1.set_ylabel("Residual")
plt.suptitle("Voigt ---   %s (Temp: %s Pressure: %s)" %(Molecule, str(int(TempValue)),str(round(P_Value,3))))
plt.tight_layout()
#plt.show()
plt.savefig("BenchmarkTest_%s_Voigt.png" %(Molecule))

plt.close()

import numpy as np
import time
import hapi
from HAPILite import CalcCrossSection
from lib.ReadComputeFunc import ReadData
import matplotlib.pyplot as plt

hapi.db_begin("TestData")

Molecule = "H2"
TempValue = 900.0
P_Value = 10.0

#TempValue = 1000.0
#P_Value = 10.0
OmegaWingValue = 100.0
OmegaRangeValue = [333.333,33333.33333]
WaveNumber = np.arange(OmegaRangeValue[0], OmegaRangeValue[1]+0.01, 0.01)
Env= {'T': TempValue, 'p': P_Value}


input("Wait here before calculation of the cross-section itself...")
#Caculate from the HAPI
nu_hapi_Voigt, abs_hapi_Voigt = hapi.absorptionCoefficient_Voigt(SourceTables=Molecule,OmegaGrid=WaveNumber,
                    HITRAN_units=False, Environment=Env,  GammaL='gamma_self', OmegaWing=OmegaWingValue, OmegaWingHW=0.0, LineShift=False)#, OmegaStep=0.01)

#
#Read the data from the location

#parse the parameters.ini which contains the information
Data = [f.split(":") for f in open("CrossSectionParams/Parameters.ini",'r+')][1:]
Values = [Item[1].split("#")[0] for Item in Data]


#Load the parameters for creating
TempStart = float(Values[0])    #Step size of the temperature
TempStop = float(Values[1])     #Step size of the temperature
TempStep = float(Values[2])     #Step size of the temperature

expP_Start = float(Values[3])   #The largest log10(pressure) in atm
expP_Stop = float(Values[4])    #The smallest log10(pressure) in atm
expP_Step = float(Values[5])    #Step size of the pressure

#Get the Index


input("Wait here for the data of the matrix....")

print("The matrix is givne")




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

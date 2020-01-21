import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from numba import jit, njit
#@jit(parallel=True,fastmath=True)

def GetWaveNumbers(StartWavelength=300, EndWavelength=30000, Resolution=100000):
        '''
        Returns the wavelengths corresponding to the resolution and in units
        of cm.
        '''
        WaveLengthValues = []

        #Converting values to
        StartWavelength = StartWavelength*1.0e-7       #nm to cm
        EndWavelength = EndWavelength*1.0e-7           #nm to cm
        WaveLengthValues = [StartWavelength]

        while WaveLengthValues[-1]<EndWavelength:
            WaveLengthValues.append(WaveLengthValues[-1]+WaveLengthValues[-1]/Resolution)
        WaveLengthValues = np.array(WaveLengthValues)
        WaveNumberRange = 1./WaveLengthValues

        WaveNumberRange = np.array(sorted(WaveNumberRange))
        return WaveLengthValues, WaveNumberRange


def SymplecticInterpolation(nu_HR,abs_HR,WaveNumberGrid):
      '''
      This function takes a cross-section at high resolution:
      nu_HR is the wavenumber in increasing order
      abs_HR is the absorption cross-section in an increasing order
      The stepsize in the WaveNumberGrid is not expected to be the equal
      '''
      #Assert the Wavenumber is strictly increasing
      IndexOrder = np.argsort(np.array(WaveNumberGrid))
      assert IndexOrder[-1]>IndexOrder[0], "The wavenumber grid should be strictly increasing."
      assert nu_HR[-1]>nu_HR[0], "The high resolution wavenumber should also be strictly increasing."

      InterpolatedValues = np.empty(len(WaveNumberGrid))

      for i in range(len(WaveNumberGrid)):
          if i+1<len(WaveNumberGrid):
              StepSize= WaveNumberGrid[i+1] - WaveNumberGrid[i]
          SelectIndex = np.abs(nu_HR-WaveNumberGrid[i])<StepSize/2.0
          InterpolatedValues[i] =  np.nanmean(abs_HR[SelectIndex])
      return InterpolatedValues

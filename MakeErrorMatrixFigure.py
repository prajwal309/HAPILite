import numpy as np
import matplotlib.pyplot as plt


TempRange = np.arange(100,901,100)
expPRange = np.arange(-5.,3.1,0.5)
PValue = 10**expPRange
#Load the data
ErrorData = np.loadtxt("ErrorValue.csv", delimiter=",")


fig, ax = plt.subplots(figsize=(16,8))
Image = plt.imshow(ErrorData, origin="lower")
plt.xlabel("Temperature (K)", fontsize=25)
plt.ylabel("Pressure (Pa)", fontsize=25)
XTicks = np.arange(0,len(PValue),1)
YTicks = np.arange(0,len(TempRange),1)
plt.xticks(XTicks, expPRange)
plt.yticks(YTicks, TempRange)
plt.colorbar(Image,fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig("ComparisonError.png")
plt.close()

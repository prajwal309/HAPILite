import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
import re

import matplotlib as mpl
mpl.rc('font',**{'sans-serif':['Helvetica'], 'size':15,'weight':'bold'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1.5})
mpl.rc('ytick',**{'major.pad':22, 'color':'k'})
mpl.rc('xtick',**{'major.pad':10,})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm','bf':'monospace:bold'})
mpl.rc('text', **{'usetex':True})
mpl.rc('text.latex',preamble=r'\usepackage{cmbright},\usepackage{relsize},'+r'\usepackage{upgreek}, \usepackage{amsmath}')
mpl.rc('contour', **{'negative_linestyle':'solid'})


def ProcessData(DataFiles):
    Resolution = 3000
    Matrix = np.zeros((Resolution,len(DataFiles)))
    WavelengthMap = np.linspace(600,28000,Resolution)

    for Row,File in enumerate(DataFiles):
        Data = open(File,'r').readlines()
        LineCenter = np.array([float(Item[3:15]) for Item in Data])
        LineIntensity = np.array([float(Item[16:26]) for Item in Data])

        WaveLength = 1./LineCenter*1e7


        Interpolator = interp1d(WaveLength, LineIntensity, fill_value=0.0, bounds_error=False)
        InterpolatedValues = Interpolator(WavelengthMap)
        Matrix[:,Row] = InterpolatedValues/(0.1*np.max(InterpolatedValues))



    Matrix[Matrix == 0] = 1e-3
    return Matrix

#Get the number of molecules
DataFiles = glob.glob("AllData/*.data")
NumMolecules = len(DataFiles)

MoleculeName = [Item.split("/")[-1][:-5] for Item in DataFiles]

#Fix the names
MoleculeLabels = []
for Name in MoleculeName:
    Integers = ["2","3","4","5","6","7","8","9"]
    TempName = Name
    for Item in Integers:
        TempName = TempName.replace(Item, "$_{\\rm "+Item+"}$")
    MoleculeLabels.append(TempName)
    print(Name, TempName)

print("MoleculeLabels::", MoleculeLabels)


Matrix = ProcessData(DataFiles)



print("The number of molecules is given by::", NumMolecules)


fig, ax = plt.subplots(figsize=(12,6.5))
ytickValues = np.arange(0.5,0.5+NumMolecules,1 )
ax.imshow(np.log(Matrix.T), origin='left', cmap='gray_r', extent=[0,3000,0,1500])
UnitY_Tick = 1500/NumMolecules
YTicks = np.arange(UnitY_Tick/2,UnitY_Tick*NumMolecules,UnitY_Tick)
ax.set_yticks(YTicks)
ax.set_yticklabels(MoleculeLabels)


SpacingX = (28000-600)/6
XTicksLabels = [int(Value) for Value in np.arange(600,28000+100+SpacingX,SpacingX)]
ax.set_xticklabels(XTicksLabels)
ax.set_xlabel("WaveLength ($\\alpha$m)")


plt.tick_params(which="both", direction="in")
plt.tight_layout()
plt.savefig("ObservationWindow.png")


PdfFigure = PdfPages("ObservationWindow.pdf")
PdfFigure.savefig(fig)
PdfFigure.close()

import sys
sys.path.append('utils')
import fileio as io
import numpy as np
import matplotlib.pyplot as plt
from plotsUtil import *

filePP = 'PPreaction.txt'

A1 = io.readMultiColFile(filePP)

fig=plt.figure()
plt.plot(A1[:,0],A1[:,1],linewidth=3,color='k',label='PelePhysics')
prettyLabels('time[s]','T[K]',14)
plotLegend()
plt.show()

import sys
sys.path.append('utils')
import fileio as io
import numpy as np
import matplotlib.pyplot as plt
from plotsUtil import *

filePP = 'PPreaction.txt'
fileCT = 'CanteraReaction_2step.txt'

A1 = io.readMultiColFile(filePP)
A2 = io.readMultiColFile(fileCT)

fig=plt.figure()
plt.plot(A1[:,0],A1[:,1],linewidth=3,color='k',label='PelePhysics')
plt.plot(A2[:,0],A2[:,1],'x',linewidth=3, color='b',label='Cantera')
prettyLabels('time[s]','T[K]',14)
plotLegend()
plt.show()


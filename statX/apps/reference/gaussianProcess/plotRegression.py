#! -*- using: utf-8 -*-

import math
import numpy as np
from matplotlib import pyplot as plt
from cProfile import label

#! the trainning data generating function
def generate( x ):
    return np.sin( 2 * math.pi * x )

#! read the training data
xTrain = np.loadtxt( './gpCoordinates.dat', dtype=np.float64, skiprows=1 )
yTrain = []
for line in open("gpHps.dat"):
    if "gaussianProcess.w" in line:
        yTrain.append( np.float64( line.split()[1] ) )
yTrain = np.array(yTrain)

#! read the prediction data
xPred     = np.loadtxt('./predictionCoord.dat', dtype=np.float64)
yPred     = np.loadtxt('./predictionMean.dat', dtype=np.float64)
sigmaPred = np.sqrt(np.loadtxt('./predictionVar.dat', dtype=np.float64))

#! define csml colors
csmlColor = {'red': '#c40000',
             'blue': '#006efe',
             'black':'#000000',
             'yellow': '#cece00',
             'darkBlue': "#004989",
             }

#! plot the regression result
plt.figure(figsize=(8,5)) 
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman'], 'size':14})

plt.plot(xPred, generate(xPred), 'k--', \
         label='Ground truth: '+r'$\sin(2\pi x)$', zorder=4)
plt.scatter(xTrain, yTrain, c=csmlColor['red'], s=50, marker='x', \
            label='Observations', zorder=3)
plt.plot(xPred, yPred, color=csmlColor['darkBlue'], ls='-', \
         label='Prediction', zorder=2)
plt.fill_between(xPred, \
                 yPred - 1.96 * sigmaPred, \
                 yPred + 1.96 * sigmaPred, \
                 alpha=.3, fc=csmlColor['darkBlue'], ec='None',  \
                 label=r'$95\%$'+' confidence interval', zorder=1)

plt.xlabel('Location ' + r'$x$')
plt.ylabel('Function value '+r'$f(x)$')
plt.legend(loc='lower left',frameon=False)
plt.savefig('GPresult.pdf')

plt.show()

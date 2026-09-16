# -*- using: ustf-8 -*-
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np


def plot(input):
    # define csml colors
    csmlColor = {'red': '#c40000',
                'blue': '#006efe',
                'black':'#000000',
                'yellow': '#cece00',
                'darkBlue': "#004989",
                }

    #! plot settings
    plt.figure(figsize=(8,5)) 
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Computer Modern Roman'], 'size':14})
    plt.xlabel(r'$x$')
    plt.ylabel(r'$p(x)$')

    #! plot predefined distribution
    mean, stdDev = 0.0, 0.0
    for line in open('input.dat'):
        if not len(line):
            continue
        elif len(line.split()) > 0 and 'mean' == line.split()[0]:
            mean = np.float64(line.split()[1])
        elif len(line.split()) > 0 and 'stdDev' == line.split()[0]:
            stdDev = np.float64(line.split()[1])
    halfFigRange = 4.0*stdDev
    normDist = stats.norm(mean, stdDev)
    x = np.linspace(mean-2.0*halfFigRange, mean+2.0*halfFigRange, 200) 
    plt.plot(x, normDist.pdf(x), color='black', ls='--', lw=2, \
             label='predefined distribution', zorder=2)
        
    #! plot the histogram of samples
    numBins = 100  
    sampleData = np.loadtxt(input[0], dtype=np.float64)
    plt.hist(sampleData, numBins, facecolor=csmlColor['blue'], density='True', alpha=0.5, \
             label='sample', zorder=1)

    #! plot the histogram
    plt.xlim(mean-halfFigRange, mean+halfFigRange)
    plt.legend(frameon=False)
    plt.savefig("histogram.pdf")    
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: plotHistogram.py  <MCMCCHAIN_FILE_NAME> \n')   
        exit()
    else:
        plot(sys.argv[1:])

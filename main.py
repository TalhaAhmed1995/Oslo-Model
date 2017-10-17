import numpy as np
from matplotlib import pyplot as plt
import random

import sites as st
import oslo as os
import heightsAnalysis as hAls
import avalanchesAnalysis as aAls

"""This is the main module, where everything is run. It imports 4 modules:
1) sites.py
2) oslo.py
1) heightsAnalysis.py
2) avalanchesAnalysis.py

All plots and calculated values as featured in the report cannot be replicated.
However, similar outputs can be produced here.
Simply clicking run here will present the results.
It is recommended that you keep the largest power variable on line 24 no larger than 7.
Also keep N <= 10,000, as larger values can take long to run.
NOTE: AFTER RUNNING, IT CAN TAKE UP TO 10 MINUTES FOR RESULTS TO APPEAR FOR SIZES L <= 128 and N <= 10,000."""

### SETTING THE HIGHEST AND LOWEST SYSTEM SIZES USING EXPONENTS FOR 2**X ###
smallestPower = 3
largestPower = 6

### SETTING THE MAIN PARAMETERS FOR THE SIMULATION ###
N = 10000 # N value used for data binning
p = 0.5 # Probability of threshold slope change
a = 1.2 # The log binning factor

### STORING THE SYSTEM OBJECTS AND ASSOCIATED TRANSIENT DURATIONS IN LISTS AND ARRAYS ###
systems = [os.Oslo(2**x) for x in range(int(smallestPower), int(largestPower)+1)]
transients = np.zeros(len(systems), dtype='int64')

def transientDurations(systems):
    """ Calculates and appends the theoretical transient durations to an array """

    for i in range(len(systems)):
        L = systems[i].getSize() # Retrieiving the size of the system
        transientTime = L**2 + L # Calculating the transient duration
        transients[i] = int(transientTime)  # Appending it to the array above

def simulateModel(systems, transients, N, p):
    """ Simulates the model by initialised, driving, relaxing and updating threshold slopes """

    for system in systems:
        system.prepareSystem() # Preparing all systems with 0 slopes

    largestTran = transients[-1] # Obtaining the transient of the largest system
    t0 = largestTran + 1000 # Going 1000 steps above to ensure steady state is reached
    tMax = t0 + int(N) # Calculating max time to run simulation

    systemHeights = [] # Storing the heights of the piles for all systems
    time = [0]

    for system in systems:
        print system.getSize()
        heights = []
        firstSite = system.getSites()[0] # Getting the first site to monitor the height
        heights.append(firstSite.getHeight())
        for i in range(1, tMax+1):
            if systems.index(system) == 0:
                time.append(i) # Each grain added accounts to a time increment
            system.driveSystem() #Adding a grain to the system
            system.relaxSystem(p) # Relaxing all sites and changing thresholds with prob p
            height = firstSite.getHeight() # Height of pile is height of first site
            heights.append(height)
        systemHeights.append(heights)

    # Converting the times and heights to numpy arrays for easier manipulation
    systemHeights = np.asarray(systemHeights, dtype='int64')
    time = np.asarray(time, dtype='int64')

    return time, systemHeights # Returning the times and heights to statistics to be performed later on

### SIMULATING THE MODEL FOR ALL SYSTEM SIZES USING THE SPECIFIED PARAMETERS ###
transientDurations(systems) # Calculating transients
time, systemHeights = simulateModel(systems=systems, transients=transients, N=N, p=p) # Obtaining the times and heights after simulating

### PERFORMING HEIGHTS ANALYSIS FOR SECTION 1.1.2 OF SCRIPT ALONG WITH PLOTS ###
heightsAn = hAls.Heights(systems, transients, time, systemHeights, N) # Creating a heights analysis object
# heightsAn.plotHeights() # Plotting heights vs time
# heightsAn.plotCrossOver() # Plotting transient durations vs L
tempTime, tempSystemHeights = heightsAn.movingAverage(w=25, plot=True) # Calculating moving average heights and times
heightsAn.dataCollapse(tempTime=tempTime, tempSystemHeights=tempSystemHeights) # Data collapsing the heights and time graphs
averageHeights, scaledHeights, a0, b, w1 = heightsAn.averageHeights(optimize=True, plot=True, plotScaled=True) # Calculating corrections to scaling parameters
STDs, scaledSTDs, a0s, Ds = heightsAn.std(optimize=True, plot=True, plotScaled=False, log=False) # Calculating STDs of average heights along with fit parameters
allHs, allProbs = heightsAn.gaussians(checkSum=False, plot=True) # Calculating height probabilities
heightsAn.collapseGaussian(allHs, allProbs, averageHeights, STDs) # Collapsing the height probabilities for all system sizes

### PERFORMING AVALANCHES ANALYSIS FOR SECTION 1.1.3 OF SCRIPT ALONG WITH PLOTS###
avAn = aAls.Avalanches(systems, transients, time, N, a) # Creating an avalanches analysis object
S, SProbs, systemCentres, systemProbs = avAn.avProb(raw=True, binned=True, plot=True, log=True) # Calculating av probs for both raw and binned
S, SProbs, systemCentres, systemProbs = avAn.avProb(raw=False, binned=True, plot=True, log=True) # Calculating only binned av probs
a, tau = avAn.getTau(systemCentres, systemProbs, plot=True, log=True) # Calculating the avalanche size exponent tau from finite-size scaling ansatz
b, D, scaledProbs, sMax = avAn.getD(systemCentres, systemProbs, tau, plot=True, plotSMax=True, log=True, log2=False) # Calculating avalanche dimension D
checkedRelation = avAn.avCollapse(systemCentres, scaledProbs, tau, D, checkValue=True, plotType='-', log=True) # Data collapsing log-binned data
kD, kTau = avAn.moments(plot=True, log=True, log2=False) # Calculating kth moments and producting plots

### PRINTING ALL THE KEY RESULTS ###
print "Transient Times:", transients
print "Maximum Time:", time[-1]
print "Steady State Heights:", averageHeights
print "Steady State STDs:", STDs
print "Average Heights a0:", a0
print "Average Heights Correction Exponent w1:", w1
print "Average Heights Correction Constant b:", b
print "Steady State STDs Scaling Relation a0:", a0s
print "Steady State STDs Scaling Relation D:", Ds
print "Values of s at Bump Peaks:", sMax
print "Avalanche Scaling Ansatz Tau:", tau
print "Avalanche Scaling Ansatz D:", D
print "Theoretical Scaling Relation for D and Tau:", checkedRelation
print "Tau Scaling Constant a:", a
print "D Scaling Constant b:", b
print "kth Moment D:", kD
print "kth Moment Tau:", kTau

plt.show()

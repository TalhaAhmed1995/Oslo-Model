import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import random

import sites as st
import oslo as os
import functions as funcs
import plotting as pt
import log_bin_CN_2016 as lb

class Avalanches:

    """A class containing all methods for the Avalanche analysis.
    This includes plotting avalanche probabilities and data collapses."""

    def __init__(self, systems, transients, time, N, a):

        self.systems = systems
        self.transients = transients
        self.time = time
        self.N = N
        self.a = a
        self.tMax = time[-1]
        self.systemSizes = np.array([system.getSize() for system in self.systems], dtype='int64')
        self.systemAvSizes = np.array([system.getAvSizes() for system in self.systems], dtype='int64')

    def avProb(self, raw=False, binned=False, checkSum=False, plot=False, plotType1='o', plotType2='o-', log=False):
        """ Calculates the probabilities of getting avalanche sizes s, P(s; L) for all systems. Also log-bins the data."""

        S = []
        SProbs = []
        systemCentres = []
        systemProbs = []

        t0 = int(self.tMax - self.N) - 1
        T = self.tMax - t0 # Setting the number of data points to use for calculating the probabilities
        lowerLimit = t0 + 1 # The lower index for the heights list to scan
        upperLimit = t0 + T # The upper index for the heights list to scan

        if raw: # Raw probs calculated only is specified
            for i in range(len(self.systems)):
                s = []
                probs = []
                avalanches = self.systemAvSizes[i] # Getting the list of avalanche sizes for the current system size
                steadyAvs = np.sort(avalanches[lowerLimit:upperLimit+1]) # Only need to scan the avalanches in the steady state
                for x in range(min(steadyAvs), max(steadyAvs)+1): # Going in ascending order of the avs list
                    s.append(x)
                    count = 0
                    for y in steadyAvs:
                        if y == x:
                            count += 1 # Counting how many times the current avalanche size occurs in the list
                        elif count >= 1:
                            break # Since the list is sorted in ascending order, if the count is >= 1 and isn't increasing further, there are no more avalanches of this size
                        else:
                            pass # If the avalanche size hasn't been found yet, keep scanning
                    prob = count / float(len(steadyAvs)) # Dividing the counted occurrences of all s values by the total number of avalanches in the steady state
                    probs.append(prob)
                s = np.asarray(s, dtype='int64') # Converting the lists to numpy arrays for easier manipulation
                probs = np.asarray(probs, dtype='float64')
                S.append(s)
                SProbs.append(probs)

            # Converting the lists to numpy arrays for easier manipulation
            S = np.asarray(S)
            SProbs = np.asarray(SProbs)

            if checkSum: # Checking to see if all the probabilities add up to 1 ie they are normalised for all system sizes
                for probs in SProbs:
                    summation = np.sum(probs)
                    print summation

        if binned: # Binned probs calculated only is specified
            for i in range(len(self.systems)):
                avalanches = self.systemAvSizes[i] # Getting the list of avalanche sizes for the current system size
                steadyAvs = avalanches[lowerLimit:upperLimit+1] # Only need to scan the avalanches in the steady state
                centres, binProbs = lb.log_bin(steadyAvs, a=self.a) # Calculating the log-binned probs using the log_bin_CN_2016.py module
                systemCentres.append(centres)
                systemProbs.append(binProbs)

        if plot: # Plotting the probs P(h; L) against avalanche size s on log-log plots (if specified)
            xLab = r"$s$"
            title = "Avalanche Size Probability vs Avalanche Size"
            if raw and not binned: # Plotting just the raw data
                yLab = r"$P_{N}(s; L)$"
                legend = ['Raw for L = ' +str(i) for i in self.systemSizes]
                plotTypes = [plotType1 for i in range(len(SProbs))]
                pt.plot(xData=S, yData=SProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)
            elif binned and not raw: # Plotting just the binned data
                yLab = r"$\tildeP_{N}(s; L)$"
                legend = ['Binned for L = ' +str(i) for i in self.systemSizes]
                plotTypes = [plotType2 for i in range(len(systemProbs))]
                pt.plot(xData=systemCentres, yData=systemProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)
            elif raw and binned: # Plotting both the raw data and binned data
                yLab = r"$P_{N}(s; L)$"
                legend = ['Raw for L = ' +str(self.systemSizes[-1]), 'Binned for L = ' +str(self.systemSizes[-1])]
                plotTypes = [plotType1, plotType2]
                pt.plot(xData=[S[-1], systemCentres[-1]], yData=[SProbs[-1], systemProbs[-1]], plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)

        return S, SProbs, systemCentres, systemProbs # Returning both the raw and binned data

    def getTau(self, systemCentres, systemProbs, scale=False, plot=False, plotType1='o', plotType2='-', log=False):
        """ Calculates the avalanche-size exponent, tau, from the finite-size scaling ansatz """

        # The probs of the largest system size vs s are optimised using a logarithmic fit function "(-tau * s) + np.log10(a)" in the linear region of the log graphs
        popt, pcov = curve_fit(funcs.func3, np.log10(systemCentres[-1][6:41]), np.log10(systemProbs[-1][6:41]))

        a, tau = popt[0], popt[1] # Getting the paramaters from the optimised fit

        if plot: # Plotting the probs P(h; L) against avalanche size s on log-log (if specified) plot for the largest system. The optimised site is also plot
            xLab = r"$s$"
            yLab = r"$\tildeP_{N}(s; L)$"
            title = "Avalanche Size Probability vs Avalanche Size (for Largest System)"
            legend = ['Binned for L = ' +str(self.systemSizes[-1]), 'Optimized fit ' r'$as^{-\tau_s}$']
            plotTypes = [plotType1, plotType2]
            pt.plot(xData=[systemCentres[-1], systemCentres[-1][6:41]], yData=[systemProbs[-1], funcs.func4(systemCentres[-1][6:41], a, tau)], plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)

        return a, tau # Returning the calculated parameter values

    def getD(self, systemCentres, systemProbs, tau, plot=False, plotSMax=False, plotType1='o-', plotType2='o', plotType3='-', log=False, log2=False):
        """ Calculates the avalanche dimension, D, from the finite-size scaling ansatz """

        sMax = [] #s values of the peaks of the bumps beyond the cutoff regions from the rescaled data
        scaledProbs = [] #s^(tau) * P(s; L)

        for i in range(len(self.systems)):
            centres = systemCentres[i]
            binProbs = systemProbs[i]
            sPower = centres ** tau # Calculating s^(tau)
            binProbs = sPower * binProbs # Calculating s^(tau) * P(s; L) (vertically collapses the data)
            index = 5
            for x in range(6, len(binProbs)-1): # Finding the s values of the peaks of the bumps beyond the cutoff regions from the rescaled data
                if binProbs[x] > binProbs[x-1]:
                    index = x
                else:
                    pass
            sMaxValue = centres[index]
            scaledProbs.append(binProbs)
            sMax.append(sMaxValue)

        # Converting the lists to numpy arrays for easier manipulation
        sMax = np.asarray(sMax)
        scaledProbs = np.asarray(scaledProbs)

        # Optimising a plot of log(sMax) vs log(L) using the equation  "(D * L) + b" (two smallest system sizes ignored to avoid corrections to scaling)
        popt, pcov = curve_fit(funcs.func5, np.log10(self.systemSizes[2:]), np.log10(sMax[2:]))
        b, D = 10**(popt[0]), popt[1]

        if plot: # Plotting the rescaled probs s^(tau) * P(s; L) vs s
            xLab = r"$s$"
            yLab = r"$s^{\tau_{s}}\tildeP_{N}(s; L)$"
            title = "Scaled Avalanche Size Probability vs Avalanche Size"
            legend = ['Scaled Binned for L = ' +str(i) for i in self.systemSizes]
            plotTypes = [plotType1 for i in range(len(scaledProbs))]
            pt.plot(xData=systemCentres, yData=scaledProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)
        if plotSMax: # Plotting sMax vs L
            xLab = r"$L$"
            yLab = r"$s_{max}$"
            title = "Avalanche Size at Bump Peak vs System Size"
            legend = ['Avalanche Size at Bump Peak', 'Optimized fit ' r'$bL^{D}$']
            plotTypes = [plotType2, plotType3]
            optimizeData = [i for i in range(self.systemSizes[0], self.systemSizes[-1]+1)]
            pt.plot(xData=[self.systemSizes, optimizeData], yData=[sMax, funcs.func6(optimizeData, b, D)], plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, log=log2)

        return b, D, scaledProbs, sMax # Returning the scaling parameters as well as the rescaled probs

    def avCollapse(self, systemCentres, scaledProbs, tau, D, checkValue=False, plotType='o-', log=False):
        """ Completely collapses the probabilities calculated using the values of tau and D calculated from the finite-size scaling ansatz """

        scaledCentres = []

        for i in range(len(self.systems)):
            centres = systemCentres[i]
            L = self.systemSizes[i]
            centres = centres / L**D # Rescaling the horizontal axis as s/L^D
            scaledCentres.append(centres)

        # Converting the list to numpy array for easier manipulation
        scaledCentres = np.asarray(scaledCentres)

        # Plotting the rescaled probs s^(tau) * P(s; L) vs rescaled avalanche sizes s/L^D for a complete data collapse
        xLab = r"$s/L^{D}$"
        yLab = r"$s^{\tau_{s}}\tildeP_{N}(s; L)$"
        title = "Scaled Avalanche Size Probability vs Scaled Avalanche Size"
        legend = ['Scaled Binned for L = ' +str(i) for i in self.systemSizes]
        plotTypes = [plotType for i in range(len(scaledProbs))]
        pt.plot(xData=scaledCentres, yData=scaledProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)

        if checkValue: # Checking if the measured values of tau and D obey the relation D*(2 - tau) = 1
            checkedRelation = D * (2 - tau)
            return checkedRelation

    def moments(self, plot=False, plotType1='o-', plotType2='o', plotType3='-', log=False, log2=False):
        """ Calculates and plots the kth moments using the avalanche sizes. Also returns estimated values of tau and D """

        kMoments = [[] for i in range(5)] # Storing the k=1,2,3,4,5 moments for all system sizes
        grads = [] # List of gradients calculated for each moment
        cs = [] # Intercepts of each moment

        t0 = int(self.tMax - self.N) - 1
        T = self.tMax - t0 # Setting the number of data points to use for calculating the kth moments
        lowerLimit = t0 + 1 # The lower index for the heights list to scan
        upperLimit = t0 + T # The upper index for the heights list to scan

        for i in range(len(self.systems)):
            avalanches = self.systemAvSizes[i] # Getting the list of avalanche sizes for the current system size
            steadyAvs = np.array(avalanches[lowerLimit:upperLimit+1], dtype='float64') # Only need to scan the avalanches in the steady state
            kMoment = []
            for k in range(1, 6): # Calculating all kth moments ranging from k=1 to 5
                values = steadyAvs ** k
                summation = np.sum(values)
                moment = summation / float(T)
                kMoments[k-1].append(moment)

        # Converting the list to numpy array for easier manipulation
        kMoments = np.asarray(kMoments)

        for moments in kMoments:
            x = [np.log10(xValue) for xValue in self.systemSizes[1:]] # Calculating log(L), ignoring the smallest system size to avoid corrections to scaling
            y = [np.log10(yValue) for yValue in moments[1:]] # Calculating log(<s^k>), again ignoring the smallest system size
            A = np.vstack([x, np.ones(len(x))]).T
            m, c = np.linalg.lstsq(A, y)[0] # Least squares regression fit to the linear set of data points
            grads.append(m) # Storing the gradients
            cs.append(c) # And the interceps

        # Converting the list to numpy array for easier manipulation
        grads = np.asarray(grads, dtype='float64')
        cs = np.asarray(cs, dtype='float64')

        ks = np.array([i+1 for i in range(len(kMoments))], dtype='int64') # Storing the k values from 1 to 5 in an array
        A = np.vstack([ks, np.ones(len(ks))]).T
        D, c = np.linalg.lstsq(A, grads)[0] # Least squares regression fit to the linear set of points of D(1 + k - tau) vs k. D is the gradient of the line.

        kValues =[]
        extrap = []
        for k in np.arange(0.1, len(kMoments)+0.01, 0.01): # Extrapolating the fit using the calculated gradient, D, and intercept, c
            y = D*k + c
            extrap.append(y)
            kValues.append(k)

        # Converting the list to numpy array for easier manipulation
        kValues = np.asarray(kValues)
        extrap = np.asarray(extrap)

        kIntercept = kValues[np.argmin(np.absolute(extrap))] # Finding the intercept of the line on the k axis
        tau = kIntercept + 1 # tau is calculated using the intercept

        if plot: # Plotting <s^k> vs L
            xLab = r"$L$"
            yLab = r"$< s^{k} >$"
            title = "kth Moment vs System Size"
            legend = ['k = ' +str(i) for i in range(1, len(kMoments)+1)]
            plotTypes = [plotType1 for i in range(len(kMoments))]
            pt.plot(xData=[self.systemSizes for i in range(len(kMoments))], yData=kMoments, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, log=log)

            # Also plotting D(1 + k - tau) vs k
            xLab = r"$k$"
            yLab = r"$D(1 + k - \tau_{s})$"
            title = "Estimated Exponent vs kth Moment"
            legend = ['Gradients', 'Extrapolated fit']
            plotTypes = [plotType2, plotType3]
            pt.plot(xData=[ks, kValues], yData=[grads, extrap], plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, log=log2)

        return D, tau # Returning the values of D and tau calculated using kth moment analysis

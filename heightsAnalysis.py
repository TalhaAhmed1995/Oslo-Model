import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import random

import sites as st
import oslo as os
import functions as funcs
import plotting as pt

class Heights:

    """A class containing all methods for the Heights analysis.
    This includes plotting heights, STDs, data collapses and probabilities."""

    def __init__(self, systems, transients, time, heights, N):

        self.systems = systems
        self.transients = transients
        self.time = time
        self.systemHeights = heights
        self.N = N
        self.tMax = time[-1]
        self.systemSizes = np.array([system.getSize() for system in self.systems], dtype='int64')

    def plotHeights(self, plotType='-', log=False):
        """ Plots the height of the pile against time for all systems """

        # Plotting the data
        xLab = r"$t$"
        yLab = r"$h(t; L)$"
        title = "Height of Pile (at Site 1) vs Time"
        legend = ['L = ' +str(i) for i in self.systemSizes]
        pt.plot(xData=self.time, yData=self.systemHeights, plotType=plotType, xLab=xLab, yLab=yLab, title=title, legend=legend, multiY=True, log=log)

    def plotCrossOver(self, plotType1='o', plotType2='-', log=False):
        """ Plots the theoretical cross-over time against L for all systems """

        # Plotting the data
        xLab = r"$L$"
        yLab = r"$t_c(L)$"
        title = "Time to Reach Steady State vs System Size"
        legend = ["Cross Over Times", r"$L^2 + L$"]
        plotTypes = [plotType1, plotType2]
        xValues = np.array([i for i in range(self.systemSizes[0], self.systemSizes[-1] + 1)], dtype='int64')
        pt.plot(xData=[self.systemSizes, xValues], yData=[self.transients, xValues**2 + xValues], plotType=plotTypes, xLab=xLab, yLab=yLab, \
        title=title, legend=legend, multiX=True, multiY=True, log=log)

    def movingAverage(self, w, plot=False, plotType='-', log=False):
        """ Smooths the heights data using a temporial time window of specified width 2w+1 """

        window = 2*w + 1 # Calculating the window width

        tempSystemHeights = []
        tempTime = []

        for i in range(len(self.systemHeights)):
            heights = self.systemHeights[i]
            tempHeights = []
            for t in range(w, len(heights)-w, 2*w):
                if i==0:
                    tempTime.append(t) # Storing the temporal time in a list
                summation = np.sum(heights[t-w:t+w+1]) # Summing the heights in the window
                tempHeight = summation/float(window) # Dividing by the window width to get an averaged height
                tempHeights.append(tempHeight) # Storing the averaged height in a list
            tempSystemHeights.append(tempHeights)

        # Converting the lists to numpy arrays to easier manipulation
        tempSystemHeights = np.asarray(tempSystemHeights, dtype='float64')
        tempTime = np.asarray(tempTime, dtype='int64')

        if plot: # Plotting the temporal heights vs temporal time
            xLab = r"$t$"
            yLab = r"$\tilde{h}(t; L)$"
            title = "Height of Pile (at Site 1) vs Time, Averaged over Time Windows, " r"$2W + 1$"
            legend = ['L = ' +str(i) for i in self.systemSizes]
            pt.plot(xData=tempTime, yData=tempSystemHeights, plotType=plotType, xLab=xLab, yLab=yLab, title=title, legend=legend, multiY=True, log=log)

        return tempTime, tempSystemHeights # Returning the averaged heights and times for all the windows

    def dataCollapse(self, tempTime, tempSystemHeights, plotType='-', log=False):
        """ Produces a data collapse of the heights against time """

        scaledSystemTimes = []
        scaledSystemHeights = []

        for i in range(len(tempSystemHeights)):
            L = self.systemSizes[i]
            tempHeights = tempSystemHeights[i]
            scaledTimes = tempTime / float(L**2) # Time is rescaled as t/L^2 as time scales quadratically with L
            scaledHeights = tempHeights / float(L) # Heights are rescaled as h/L as heights scale linearly with L
            scaledSystemTimes.append(scaledTimes)
            scaledSystemHeights.append(scaledHeights)

        # Converting the lists to numpy arrays to easier manipulation
        scaledSystemTimes = np.asarray(scaledSystemTimes, dtype='float64')
        scaledSystemHeights = np.asarray(scaledSystemHeights, dtype='float64')

        # Plotting the rescaled heights vs rescaled time
        xLab = r"$t/L^2$"
        yLab = r"$\tilde{h}(t; L)/L$"
        title = "(Height of Pile)/" r"$L$" " vs Time/" r"$L^2$" " (DATA COLLAPSE)"
        legend = ['L = ' +str(i) for i in self.systemSizes]
        plotTypes = [plotType for i in range(len(scaledSystemHeights))]
        pt.plot(xData=scaledSystemTimes, yData=scaledSystemHeights, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, \
        multiY=True, loc=4, log=log)

    def averageHeights(self, optimize=False, plot=False, plotScaled=False, plotType1='o', plotType2='-', log=False):
        """ Calculates the average height of the pile once in the steady state for all systems. Also computes the corrections to scaling parameters"""

        averageHeights = []
        scaledHeights = []

        t0 = int(self.tMax - self.N)
        T = self.tMax - t0 # Setting the number of data points to use for calculating the average heights
        lowerLimit = t0 + 1 # The lower index for the heights list to scan
        upperLimit = t0 + T # The upper index for the heights list to scan

        for i in range(len(self.systemHeights)):
            heights = self.systemHeights[i]
            L = self.systemSizes[i]
            steady = heights[lowerLimit:upperLimit+1] # Getting the steady state heights
            averageHeight = np.sum(steady)/float(T) # Summing the steady state heights and dividing by T
            scaledHeight = averageHeight/float(L) # Scaling the average heights as <h>/L to see corrections to scaling
            averageHeights.append(averageHeight)
            scaledHeights.append(scaledHeight)

        # Converting the lists to numpy arrays to easier manipulation
        averageHeights = np.asarray(averageHeights, dtype='float64')
        scaledHeights = np.asarray(scaledHeights, dtype='float64')

        if optimize: # The scaled heights are optimised using a fit function "a0 - b*(L**(-w1))"
            popt, pcov = curve_fit(funcs.func1, self.systemSizes, scaledHeights)
            a0, b, w1 = popt[0], popt[1], popt[2] # Getting the paramaters for the corrections to scaling for the average heights

        if plot: # Plotting the average heights vs L
            xLab = r"$L$"
            yLab = r"$<h(t; L)>_t$"
            title = "Average Height of Pile (at Site 1) in Steady State vs System Size"
            legend = "Steady State Heights"
            pt.plot(xData=self.systemSizes, yData=averageHeights, plotType='o-', xLab=xLab, yLab=yLab, title=title, legend=legend, log=log)

        if plotScaled: # Plotting the rescaled heights <h>/L vs L to see signatures of corrections to scaling
            xLab = r"$L$"
            yLab = r"$<h(t; L)>_{t}/L$"
            title = "(Average Height of Pile in Steady State)/" r"$L$" " vs System Size"
            legend = ["Scaled Steady State Heights", "Optimized Fit Function " r"$a_0 - bL^{-w_1}$"]
            plotTypes = [plotType1, plotType2]
            optimizeData = [i for i in range(self.systemSizes[0], self.systemSizes[-1]+1)]
            pt.plot(xData=[self.systemSizes, optimizeData], yData=[scaledHeights, funcs.func1(optimizeData, a0, b, w1)], plotType=plotTypes, \
            xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, log=log)

        if optimize:
            return averageHeights, scaledHeights, a0, b, w1 # If optimisation was used, the corrections to scaling parameters are also returned
        else:
            return averageHeights, scaledHeights # Otherwise just the average and scaled heights are returned


    def std(self, optimize=False, plot=False, plotScaled=False, plotType1='o', plotType2='-', log=False):
        """ Calculates the STD of the height of the pile once in the steady state for all systems. Also produces an optimised fit through the points"""

        STDs = []

        t0 = int(self.tMax - self.N)
        T = self.tMax - t0 # Setting the number of data points to use for calculating the stds
        lowerLimit = t0 + 1 # The lower index for the heights list to scan
        upperLimit = t0 + T # The upper index for the heights list to scan

        for i in range(len(self.systemHeights)):
            heights = self.systemHeights[i]
            steady = heights[lowerLimit:upperLimit+1]
            firstTerm = np.sum(steady**2)/float(T) # Calculating <h^2>
            secondTerm = (np.sum(steady)/float(T))**2 # Calculating <h>^2
            STD = np.sqrt(firstTerm - secondTerm) # Substracting the two terms and taking the sqrt
            STDs.append(STD)

        # Converting the list to numpy array to easier manipulation
        STDs = np.asarray(STDs, dtype='float64')

        if optimize:# The STDs are optimised using a fit function "a0 * (L**D)" (two smalllest system sizes ignored to avoid corrections to scaling)
            popt, pcov = curve_fit(funcs.func2, self.systemSizes[2:], STDs[2:])
            a0, D = popt[0], popt[1]  # Getting the paramaters from the optimised fit
            scaledSTDs = np.array(STDs / (self.systemSizes**float(D)), dtype='float64') # Scaling the STDs by dividing by L^D

        if plot: # Plotting STDs vs L
            xLab = r"$L$"
            yLab = r"$\sigma_{h}(L)$"
            title = "Average Height STD of Pile in Steady State vs System Size"
            legend = "Average Height STDs"
            if optimize: # Also plotting the optimised fit curve through the numerical data point
                legend = ["Average Height STDs", "Optimized Fit Function " r"$a_0L^D$"]
                plotTypes = [plotType1, plotType2]
                optimizeData = [i for i in range(self.systemSizes[0], self.systemSizes[-1]+1)]
                pt.plot(xData=[self.systemSizes, optimizeData], yData=[STDs, funcs.func2(optimizeData, a0, D)], plotType=plotTypes, xLab=xLab, yLab=yLab, \
                title=title, legend=legend, multiX=True, multiY=True, log=log)
            else:
                pt.plot(xData=self.systemSizes, yData=STDs, plotType='o-', xLab=xLab, yLab=yLab, title=title, legend=legend, log=log)

        if plotScaled: # Plotting scaled STDs vs L
            xLab = r"$L$"
            yLab = r"$\sigma_{h}(L)/L^{D}$"
            title = "(Average Height STD)/" r"$L^{D}$" " vs System Size"
            legend = "Scaled Average Height STDs"
            pt.plot(xData=self.systemSizes, yData=scaledSTDs, plotType='o-', xLab=xLab, yLab=yLab, title=title, legend=legend, log=log)

        if optimize:
            return STDs, scaledSTDs, a0, D # Also returning fit parameters
        else:
            return STDs # Only returning calculated STDs

    def gaussians(self, checkSum=False, plot=False, plotType='o-', log=False):
        """ Calculates the probabilities of getting heights h, P(h; L) for all systems """

        allHs = []
        allProbs = []

        t0 = int(self.tMax - self.N)
        T = self.tMax - t0 # Setting the number of data points to use for calculating the probabilities
        lowerLimit = t0 + 1 # The lower index for the heights list to scan
        upperLimit = t0 + T # The upper index for the heights list to scan

        for i in range(len(self.systemHeights)):
            h = []
            probs = []
            heights = self.systemHeights[i]
            steady = np.sort(heights[lowerLimit:upperLimit+1]) # Only need to scan the heights in the steady state
            for x in range(min(steady), max(steady)+1): # Going in ascending order of the heights list
                h.append(x)
                count = 0
                for y in steady:
                    if y == x:
                        count += 1 # Counting how many times the current height value occurs in the list
                    elif count >= 1:
                        break # Since the list is sorted in ascending order, if the count is >= 1 and isn't increasing further, there are no more heights of this value
                    else:
                        pass # If the height value hasn't been found yet, keep scanning
                prob = count / float(len(steady)) # Dividing the counted occurrences of all h values by the total number of heights in the steady state
                probs.append(prob)
            h = np.asarray(h, dtype='int64') # Converting the lists to numpy arrays for easier manipulation
            probs = np.asarray(probs, dtype='float64')
            allHs.append(h)
            allProbs.append(probs)

        # Converting the lists to numpy arrays for easier manipulation
        allHs = np.asarray(allHs)
        allProbs = np.asarray(allProbs)

        if checkSum: # Checking to see if all the probabilities add up to 1 ie they are normalised for all system sizes
            for probs in allProbs:
                summation = np.sum(probs)
                print summation

        if plot: # Plotting the probabilities against h for all systems
            xLab = r"$h$"
            yLab = r"$P(h; L)$"
            title = "Height Probability vs Height Value"
            legend = ['L = ' +str(i) for i in self.systemSizes]
            plotTypes = [plotType for value in range(len(allHs))]
            pt.plot(xData=allHs, yData=allProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)

        return allHs, allProbs # Returning the probs with the associated h values

    def collapseGaussian(self, Hs, probs, averageHeights, STDs, plotType='o-', log=False):
        """ Collapses the probabilities calculated into a centralised Gaussian """

        scaledHs = np.array((Hs - averageHeights)/STDs) # Scaling the h values as (h - <h>) / sigma
        scaledProbs = np.array(probs * STDs) # Scaling the probs as sigma * Probs

        # Plotting the rescaled probabilities vs the rescaled h values
        xLab = r"$(h - <h>)/\sigma_{h}$"
        yLab = r"$\sigma_{h}P(h; L)$"
        title = "Scaled Height Probability vs Scaled Height (DATA COLLAPSE)"
        legend = ['L = ' +str(i) for i in self.systemSizes]
        plotTypes = [plotType for value in range(len(scaledHs))]
        pt.plot(xData=scaledHs, yData=scaledProbs, plotType=plotTypes, xLab=xLab, yLab=yLab, title=title, legend=legend, multiX=True, multiY=True, loc=1, log=log)

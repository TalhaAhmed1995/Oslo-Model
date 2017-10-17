import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import random

"""This module is strictly used for plotting various graphs."""

def plot(xData, yData, plotType, xLab, yLab, title, legend, multiX=False, multiY=False, loc=2, log=False):
    """ Plots the input x and y data. Customised titles and labels can also be input, as well as choosing the type of scale """

    plt.figure(figsize=(12, 10))
    if multiX == True and multiY == False: # If multiple x data points need to be plotted on the same y points
        for i in range(len(xData)):
            plt.plot(xData[i], yData, plotType, label=legend[i])
    elif multiY == True and multiX == False: # If multiple y data points need to be plotted on the same x points
        for i in range(len(yData)):
            plt.plot(xData, yData[i], plotType, label=legend[i])
    elif multiX == True and multiY == True: # If multiple sets of data need to be plotted
        for i in range(len(xData)):
            plt.plot(xData[i], yData[i], plotType[i], label=legend[i])
    else:
        plt.plot(xData, yData, plotType, label=legend)
    plt.xlabel(xLab, fontsize=20)
    plt.ylabel(yLab, fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(title, fontsize=18)
    plt.legend(loc=loc, prop={'size':20})
    if log: # The axes are scaled logarithmically if specific
        plt.axes().set_xscale('log')
        plt.axes().set_yscale('log')
    plt.grid()

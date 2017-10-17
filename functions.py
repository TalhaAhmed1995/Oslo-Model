import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import random

"""This module is strictly used for the functions required to produce optimised fits and calculating parameters."""

def func1(L, a0, b, w1):
    """Calculates function for corrections to scaling for average heights"""

    return a0 - b*(L**(-w1))

def func2(L, a0, D):
    """Calculates function for optimised fit of STD"""

    return a0 * (L**D)

def func3(s, a, tau):
    """ Calculates logarithmic form of the avalanche-size exponent scaling ansatz """

    return (-tau * s) + np.log10(a)

def func4(s, a, tau):
    """ Calculates the exponent form of the avalanche-size exponent scaling ansatz """

    return a * (s**(-tau))

def func5(L, b, D):
    """ Calculates the logarithmic form of the avalanche dimension scaling ansatz """

    return (D * L) + b

def func6(L, b, D):
    """ Calculates the exponent form of the avalanche dimension scaling ansatz """

    return b * (L**D)

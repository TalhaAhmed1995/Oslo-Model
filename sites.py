import numpy as np
from matplotlib import pyplot as plt

class Site:

    """A class containing all methods for each site location of a system.
    This includes returning and updating the height, location, slope and threshold slope."""

    def __init__(self, location, threshold):

        self.location = location
        self.height = 0
        self.threshold = threshold
        self.slope = 0

    def addGrain(self):
        """ Increases the height by 1 unit """

        self.height += 1

    def removeGrain(self):
        """ Decreases the height by 1 unit """

        self.height -= 1

    def getHeight(self):
        """ Returns the current height """

        return self.height

    def setLocation(self, location):
        """ Sets the location of the site if it needs to be changed """

        self.location = location

    def getLocation(self):
        """ Returns the current location """

        return self.location

    def setThreshold(self, threshold):
        """ Sets a new threshold slope value """

        self.threshold = threshold

    def getThreshold(self):
        """ Returns the current threshold slope """

        return self.threshold

    def changeSlope(self, slope):
        """ Increases/decreases the slope by a specified value """

        self.slope += slope

    def getSlope(self):
        """ Returns the current slope """

        return self.slope

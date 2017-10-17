import numpy as np
from matplotlib import pyplot as plt
import random

import sites as st

class Oslo:

    """A class containing all methods for the Oslo model algorithm.
    This includes initialising, driving, relaxing and updating the threshold slopes of the system."""

    def __init__(self, size):

        self.size = size
        self.sites = []
        self.avSizes = []

    def getSize(self):
        """ Returns the size, L, of the system """

        return self.size

    def getAvSizes(self):
        """ Returns a list of the avalanche sizes of the system """

        return self.avSizes

    def getSites(self):
        """ Returns a list of all the site objects in the system """

        return self.sites

    def prepareSystem(self):
        """ Initialises the system with all sites set to slope zero and randomly chosen threshold slopes """

        del self.sites[:] # Deleting any previous sites in the system
        for i in range(self.size): # Creating new sites for range of system size
            location = i + 1 # Setting the location of the site
            threshold = random.randint(1, 2) # Randomly choosing the threshold slope
            site = st.Site(location, threshold) # Creating new site object
            self.sites.append(site)

    def randomThreshold(self, p):
        """ Randomly sets the threshold slope of a site. The specified value of p gives the prob of getting \
        a threshold slope of 1, while 1 - p is the prob of threshold 2 """

        randVariate = random.random() # Uniform random number between 0 and 1
        if randVariate <= p: # IF value is <= p, threshold is set to 1
            newThreshold = 1
        elif randVariate > p: # Otherwise it's set to 2
            newThreshold = 2

        return newThreshold # Returning the newly selected threshold slope value

    def driveSystem(self):
        """ Drives the system by adding a grain to site location 1 """

        if len(self.sites) == 0: # Checking if the system has been prepared first. If not, a warning is returned.
            raise Exception("Your system has not been initialised with sites. Please define the size of the system > 0.")
        else:
            firstSite = self.sites[0]  # Getting the first site
            firstSite.addGrain() # Increasing its height by 1
            firstSite.changeSlope(1) # Increasing its slope by 1

    def relaxSystem(self, p):
        """ Relaxes the entire system by checking the slope against the threshold slope for each site.\
        Each site is checked multiple times to ensure that all sites are below their threshold slopes """

        avalancheSize = 0 # Counter for the number of topplings that have occured

        allRelaxed = False
        for site in self.sites:
            if site.getSlope() > site.getThreshold():
                allRelaxed = False # If even a single site has slope > threshold, the system is deemed not relaxed
                break
            else:
                allRelaxed = True

        while allRelaxed is False: # Looping over all sites multiple time until slope < threshold is obeyed
            for site in self.sites:
                if site.getLocation() < self.size:
                    nextSite = self.sites[site.getLocation()] # Getting neighbour i + 1 for all sites from 1 to L - 1
                if site.getLocation() > 1:
                    prevSite = self.sites[site.getLocation()-2] # Getting neighbour i - 1 for all sites from 2 to L

                if site.getSlope() > site.getThreshold(): # Checking slope against threshold
                    site.removeGrain() # Grain is toppled if so
                    avalancheSize += 1 # Toppling equates to an avalanche of size 1, so avalanche counter is increased
                    if site.getLocation() < self.size:
                        site.changeSlope(-2)
                        nextSite.addGrain() # Grain moved to next site if i < L
                        nextSite.changeSlope(1)
                        site.setThreshold(self.randomThreshold(p)) # Randomly selecting a new threshold for the relaxed site with specified prob p
                        if site.getLocation() > 1:
                            prevSite.changeSlope(1)
                    elif site.getLocation() == self.size:
                        site.changeSlope(-1) # If i = L, the grain just leaves the system
                        prevSite.changeSlope(1)
                        site.setThreshold(self.randomThreshold(p)) # Randomly selecting a new threshold for the relaxed site with specified prob p

            for site in self.sites: # Checking all sites against their thresholds again to see if anymore relaxing is needed
                if site.getSlope() > site.getThreshold():
                    allRelaxed = False
                    break
                else:
                    allRelaxed = True

        self.avSizes.append(avalancheSize) # Final count for the avalanche size (after relaxing has finished) is stored in a list

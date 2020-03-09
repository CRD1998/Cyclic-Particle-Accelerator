import log
from abc import ABC, abstractmethod
import numpy as np
import math
import copy
from statistics import mean
import scipy.constants as const

class Bunch(ABC):
    """
    This abstract base class provides the skeleton for child classes that will generate a bunch of
    charged particle objects. 

    Concrete Methods:
    -----------------
        1) assignPositions - returns an array of 3darrays, sampled from a normal
        distribution with a mean of [0,0,0].
        2) distributeEnergies - returns an array of kinetic energies (in eV), sampled from a
        normal distribution with the mean (AverageKinetic) passed into the child classes.
        3) KineticEnergy - takes an optional arguement "total" (defaults to False). Returns the 
        average kinetic energy for a particle in a bunch (eV). If total is True then it returns
        the kinetic energy of the entire bunch.
        4) momentum - takes an optional arguement "total" (defaults to False). Returns the average 
        three-momentum for a particle in a bunch (kg m/s). If total is True then it returns the 
        momentum of the entire bunch.
        5) spread - returns the standard deviation of the x-positions, y-positions and z-positions of 
        all the particles in the bunch as 3D numpy array.

    Abstract Methods:
    -----------------
        1) assignVelocites - every bunch class must implement a method that returns calculates
        a linear speed from every kinetic energy that is returned from the distributeEnergies
        method. This method is abstract because it depends on the chosen particle's mass.
        2)  createBunch - every bunch class must implement a method that actually generates
        the final list of Charged Particle objects.
    """

    conversion = const.physical_constants['electron volt'][0] # eV => joules

    @abstractmethod
    def __init__(self, AverageKinetic, particleNumber=3):
        self.average_kinetic_energy= float(AverageKinetic)
        self.bunch_number = int(particleNumber)
        
        self.bunch = self.createBunch()

        super(Bunch, self).__init__()

    @abstractmethod
    def createBunch(self):
        pass

    @abstractmethod
    def assignVelocities(self):
        pass

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError('Cannot add bunches of different types together.')
        new_bunch = copy.deepcopy(self) # make a copy of the current instance 
        new_bunch.bunch += other.bunch # add the bunches of the two parsed Bunch objects together
        new_bunch.bunch_number += other.bunch_number
        new_bunch.bunchName = 'combined ' + other.bunchName
        delattr(new_bunch, 'average_kinetic_energy')
        return new_bunch # return a Bunch object made up of the two parsed bunches

    def assignPositions(self):
        mu = 0.
        sigma = 0.01*0.001 # 1% of 1cm
        positions = np.random.normal(mu, sigma, (self.bunch_number,3))
        for i in positions:
            i[2] = 0. # set all z values to zero
        return positions

    def distributeEnergies(self):
        mu = self.average_kinetic_energy
        sigma = 0.01 * mu
        return np.random.normal(mu, sigma, self.bunch_number)
    
    def averagePosition(self):
        return np.array(np.mean([i.position for i in self.bunch],axis=0),dtype=float)

    def averageVelocity(self):
        return np.array(np.mean([i.velocity for i in self.bunch],axis=0),dtype=float)

    def KineticEnergy(self, total=False):
        if total == False:
            return mean([i.KineticEnergy() for i in self.bunch])/self.conversion
        return sum([i.KineticEnergy() for i in self.bunch])/self.conversion

    def momentum(self, total=False):
        if total == False:
            return np.array(np.mean([i.momentum() for i in self.bunch],axis=0),dtype=float)
        return np.array(np.sum([i.momentum() for i in self.bunch],axis=0),dtype=float)

    def spread(self):
        return np.array(np.std([i.position for i in self.bunch],axis=0),dtype=float)
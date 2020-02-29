import log
from abc import ABC, abstractmethod
import numpy as np
import math
import scipy.constants as const

class Bunch(ABC):
    """
    This abstract base class provides the skeleton for child classes that will generate a bunch of
    charged particle objects. 

    Concrete Methods:
        1) assignPositions - returns an array of 3darrays, sampled from a normal
        distribution with a mean of [0,0,0].
        2) distributeEnergies - returns an array of kinetic energies (in eV), sampled from a
        normal distribution with the mean (AverageKinetic) passed into the child classes.
    Abstract Methods:
        3) assignVelocites - every bunch class must implement a method that returns calculates
        a linear speed from every kinetic energy that is returned from the distributeEnergies
        method. This method is abstract because it depends on the chosen particle's mass.
        4)  createBunch - every bunch class must implement a method that actually generates
        the final list of Charged Particle objects.
    """

    conversion = const.physical_constants['electron volt'][0] # eV => joules

    @abstractmethod
    def __init__(self, AverageKinetic, particleNumber=3, positionSigma=0.01, energySigma=0.01):
        self.average_kinetic_energy= AverageKinetic
        self.bunch_number = particleNumber
        self.positionSigma = positionSigma
        self.energySigma = energySigma
        
        self.bunch = self.createBunch()

        super(Bunch, self).__init__()

    @abstractmethod
    def createBunch(self):
        pass

    @abstractmethod
    def assignVelocities(self):
        pass

    def assignPositions(self):
        mu, sigma = 0., self.positionSigma
        return np.random.normal(mu, sigma, (self.bunch_number,3))

    def distributeEnergies(self):
        mu, sigma = self.average_kinetic_energy, self.energySigma
        return np.random.normal(mu, sigma, self.bunch_number)
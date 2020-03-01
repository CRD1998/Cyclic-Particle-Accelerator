import log
from abc import ABC, abstractmethod
import numpy as np
import math
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
    def __init__(self, AverageKinetic, particleNumber=3, positionSigma=0.01, energySigma=0.01):
        self.average_kinetic_energy= float(AverageKinetic)
        self.bunch_number = int(particleNumber)
        self.positionSigma = float(positionSigma)
        self.energySigma = float(energySigma)
        
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
import log
import scipy.constants as const
import math
from statistics import mean
from Bunch import Bunch
from ChargedParticle import ChargedParticle

class ProtonBunch(Bunch):
    """
    Class to generate a bunch of protons as Charged Particle objects, the list can be accessed 
    through the bunch attribute (derived from Bunch abstract base class). 

    Attributes:
    -----------
    AverageKinetic: float/int
        The average kinetic energy of a particle in the bunch, 
        measured in eV.
    
    particleNumber: int, optional 
        Number of particles in the bunch (defaults to 3).
    
    positionSigma: float/int, optional
        The S.D of the particles' normal position distribution 
        about the origin in metres (defaults to 0.01 m).

    Methods:
    --------
    createBunch
        Returns a list of Charged Particle objects representing a bunch of protons.
    
    assignVelocities
        Returns a list of the initial 3D velocity vectors for all the particles in the bunch in m/s.
    """

    def __init__(self, AverageKinetic, particleNumber=3, positionSigma=0.01):
        self.particleName = 'proton'
        self.particleMass = const.m_p
        self.particleCharge = const.e
        self.bunchName = self.particleName + '-bunch'
        super().__init__(AverageKinetic=AverageKinetic, particleNumber=particleNumber, positionSigma=positionSigma)

    def createBunch(self):
        """
        This method returns a list of charged particle objects to represent a bunch of protons. This list
        is then assigned to bunch attribute of the ProtonBunch instance.
        """

        positions = self.assignPositions() # get the positions for all the particles in the bunch
        velocities = self.assignVelocities() # get the initial velocities of particles 
        proton_bunch= []
        i = 0
        log.logger.info('generating bunch')
        for (j,k) in zip(positions,velocities):
            i += 1
            proton_bunch.append(ChargedParticle(self.particleName+'-'+str(i), self.particleMass, 
            self.particleCharge, j, k))
        log.logger.info('bunch generated')
        return proton_bunch

    def assignVelocities(self):
        """
        Using the sampled kinetic energies in the distributeEnergies method, this method will convert them
        to linear speeds and return a list of velocity vectors where the linear speeds make up the x-component
        of these velocity vectors.
        """

        energies = self.distributeEnergies() # get the kinetic energies
        gamma = lambda x: 1 + (x*self.conversion)/(self.particleMass*const.c*const.c) # function that returns gamma, given kinetic energy
        speed = lambda x: const.c*math.sqrt(1-1/(x*x)) # function that returns a linear speed, given gamma
        return [[speed(gamma(i)),0,0] for i in energies]
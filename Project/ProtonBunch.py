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
        1) AverageKinetic - the average kinetic energy of a particle in the bunch, 
        measured in eV.
    
    Optional Attributes:
    --------------------
        2) particleNumber - number of particles in the bunch, defaults to 3.
        3) positionSigma - the S.D of the particles' normal position distribution 
        about the origin, defaults to 0.01 .
        4) energySigma - the S.D of the particles' normal kinetic distribution about 
        AverageKinetic, defaults to 0.01 .
    """
    def __init__(self, AverageKinetic, particleNumber=3, positionSigma=0.01, energySigma=0.01):
        self.particleName = 'proton'
        self.particleMass = const.m_p
        self.particleCharge = const.e
        super().__init__(AverageKinetic=AverageKinetic, particleNumber=particleNumber, 
        positionSigma=positionSigma, energySigma=energySigma)

    def createBunch(self):
        positions = self.assignPositions()
        velocities = self.assignVelocities()
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
        energies = self.distributeEnergies() # get distribution of kineitc energies
        gamma = lambda x: 1 + (x*self.conversion)/(self.particleMass*const.c*const.c) # function that returns gamma, given kinetic energy
        speed = lambda x: const.c*math.sqrt(1-1/(x*x)) # function that returns linear speed, given gamma
        return [[speed(gamma(i)),0,0] for i in energies]
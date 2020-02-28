import log
from ChargedParticle import ChargedParticle
import math
import numpy as np
import scipy.constants as const

class ProtonBunch:
    """
    This class will generate a bunch of protons, as Charged Particle objects, the list can be accessed 
    through the bunch attribute.

    particle number - number of particles in the desired bunch.
    AverageKinetic - the average kinetic energy of a particle in the bunch, measured in eV.
    sigma - the S.D of the particles' kinetic energy distribution.
    """
    def __init__(self, AverageKinetic, particleNumber=10, sigma=0.01):
        self.bunch_number = particleNumber
        self.average_kinetic_energy= AverageKinetic
        self.sigma = sigma
        self.bunch = self.createBunch()
    
    def createBunch(self):
        positions = self.assignPositions()
        velocities = self.assignVelocities()
        proton_bunch= []
        i = 0
        log.logger.info('generating bunch')
        for (j,k) in zip(positions,velocities):
            i += 1
            proton_bunch.append(ChargedParticle('proton-'+str(i), const.m_p, const.e, j, k))
        log.logger.info('bunch generated')
        return proton_bunch

    def assignPositions(self):
        mu, sigma = 0., 0.01
        return np.random.normal(mu, sigma, (self.bunch_number,3))

    def assignVelocities(self):
        conversion = const.physical_constants['electron volt'][0] # eV => joules
        energies = self.distributeEnergies() # get distribution of kineitc energies
        gamma = lambda x: 1 + (x*conversion)/(const.m_p*const.c*const.c) # function that returns gamma, given kinetic energy
        speed = lambda x: const.c*math.sqrt(1-1/(x*x)) # function that returns linear speed given gamma
        return [[speed(gamma(i)),0,0] for i in energies]

    def distributeEnergies(self):
        mu, sigma = self.average_kinetic_energy, self.sigma
        return np.random.normal(mu, sigma, self.bunch_number)
            
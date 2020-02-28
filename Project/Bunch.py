import log
from ChargedParticle import ChargedParticle
import math
import numpy as np
import scipy.constants as const

class Bunch:
    """
    This class will generate a bunch of protons, as Charged Particle objects, the list can be accessed 
    through the bunch attribute.

    Arguements:
        1) particleType - a charged particle object which will determine what type of particle is in the bunch. Only the name, mass and charge are used in this class, other attributes are ignored.
        2) AverageKinetic - the average kinetic energy of a particle in the bunch, measured in eV.
    
    Optional Arguements:
        3) particle number - number of particles in the bunch, defaults to 3.
        4) positionSigma - the S.D of the particles' normal position distribution about the origin, defaults to 0.01 .
        5) energySigma - the S.D of the particles' normal kinetic distribution about AverageKinetic, defaults to 0.01 .
    """
    def __init__(self, particleType, AverageKinetic, particleNumber=10, positionSigma=0.01, 
                 energySigma=0.01):
        self.particleName = particleType.name
        self.particleMass = particleType.mass
        self.particleCharge = particleType.charge
        self.average_kinetic_energy= AverageKinetic
        self.bunch_number = particleNumber
        self.positionSigma = positionSigma
        self.energySigma = energySigma
        
        self.bunch = self.createBunch()
    
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

    def assignPositions(self):
        mu, sigma = 0., self.positionSigma
        return np.random.normal(mu, sigma, (self.bunch_number,3))

    def assignVelocities(self):
        conversion = const.physical_constants['electron volt'][0] # eV => joules
        energies = self.distributeEnergies() # get distribution of kineitc energies
        gamma = lambda x: 1 + (x*conversion)/(self.particleMass*const.c*const.c) # function that returns gamma, given kinetic energy
        speed = lambda x: const.c*math.sqrt(1-1/(x*x)) # function that returns linear speed, given gamma
        return [[speed(gamma(i)),0,0] for i in energies]

    def distributeEnergies(self):
        mu, sigma = self.average_kinetic_energy, self.energySigma
        return np.random.normal(mu, sigma, self.bunch_number)
            
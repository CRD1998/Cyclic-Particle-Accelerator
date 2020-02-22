import log
import numpy as np
import math
import copy
import scipy.constants as const

class Particle:
    """
    Class to model a particle. 
    It will make use of numpy arrays to store the position velocity etc. 
    Working directly from past exercises... 

    rest mass in kg 
    position and velocity in m 
    """

    def __init__(self,  Name='Point', Mass=1.0, Position=np.array([0,0,0], dtype=float), Velocity=np.array([0,0,0], dtype=float), Acceleration=np.array([0,0,0], dtype=float)):
        self.name = Name
        self.position = np.array(Position,dtype=float)
        self.velocity = np.array(Velocity,dtype=float)
        self.acceleration = np.array(Acceleration,dtype=float)
        self.mass = Mass

        if self.magnitude(self.velocity) >= const.c:
            log.logger.error('%s cannot have a speed greater than the speed of light' % self.name)
            raise ValueError('Velocity cannot be greater than the speed of light.')
        else:
            log.logger.info('%s has been generated' % self.name)

    def __repr__(self):
        return 'Particle: {0}, Mass: {1:12.3e}, Position: {2}, Velocity: {3}, Acceleration: {4}'.format(self.name,self.mass,self.position, self.velocity,self.acceleration)

    def gamma(self):
        """
        Returns the Lorentz factor for any given Particle object.
        """
        return 1/(math.sqrt(1-(self.magnitude(self.velocity)*self.magnitude(self.velocity))/(const.c*const.c)))

    def magnitude(self, vector):
        return np.linalg.norm(vector)

    def KineticEnergy(self):
        return 0.5*self.mass*self.magnitude(self.velocity)*self.magnitude(self.velocity)
  
    def momentum(self):
        return self.mass*np.array(self.velocity,dtype=float)

    def euler(self, deltaT):
        self.position +=  self.velocity*deltaT
        self.velocity +=  self.acceleration*deltaT

    def eulerCromer(self, deltaT):
        self.velocity +=  self.acceleration*deltaT
        self.position +=  self.velocity*deltaT

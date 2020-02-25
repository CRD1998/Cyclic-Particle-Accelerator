import log
import numpy as np
import scipy.constants as const
import math
from System import system

class EMField:

    def __init__(self, ElectricField = np.array([0,0,0], dtype=float), 
                MagneticField = np.array([0,0,0], dtype=float)):
        self.electric = np.array(ElectricField,dtype=float)
        self.magnetic = np.array(MagneticField,dtype=float)
        log.logger.info('electromagentic field generated')

    def __repr__(self):
        return 'EM Field: \n Electric Field = {0}  \n Magnetic Field = {1}'.format(self.electric,self.magnetic)

    def electricMag(self):
        return np.linalg.norm(self.electric)
    
    def magneticMag(self):
        return np.linalg.norm(self.magneticMag)
    
    def getAcceleration(self,particle):
        lorentz = np.array([0,0,0],dtype=float)
        if particle.position[0] > -0.98:
            lorentz += -1*self.electric
        if particle.position[0] < 0.098:
            lorentz += self.electric
        lorentz += np.cross(particle.velocity, self.magnetic)
        lorentz *= (particle.charge / particle.mass)
        particle.acceleration = lorentz

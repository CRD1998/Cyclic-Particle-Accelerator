import log
import numpy as np
import scipy.constants as const
import math
from System import system

class EMField:
    """
    Class to represent an  electromagnetic field.

    Parse in the electric field at its maximum and it will be converted to a time-varying
    electric field of the form Acos(wt+phi) 

    Electric field in N/C or V/m
    Magnetic field in T
    """

    def __init__(self, ElectricField = np.array([0,0,0], dtype=float), 
                MagneticField = np.array([0,0,0], dtype=float)):
        self.magnetic = np.array(MagneticField,dtype=float)
        self.electric = np.array(ElectricField,dtype=float)
        log.logger.info('electromagentic field generated')

    def __repr__(self):
        return 'EM Field: \n Electric Field = {0}  \n Magnetic Field = {1}'.format(self.electric,self.magnetic)

    def electricMag(self):
        return np.linalg.norm(self.electric)
    
    def magneticMag(self):
        return np.linalg.norm(self.magnetic)
    
    def frequency(self,particle):
        return abs(particle.charge)*self.magneticMag()/particle.mass

    def ImplementElectricField(self,particle,time):
        if particle.position[0] > -0.98 and particle.position[0] < 0.98:
            return [i*math.cos(self.frequency(particle)*time) for i in self.electric]
        return [0,0,0]

    def getAcceleration(self,particle,time):
        lorentz = np.array(self.ImplementElectricField(particle,time),dtype=float)
        lorentz += np.cross(particle.velocity, self.magnetic)
        lorentz *= particle.charge/ (particle.mass*particle.gamma())
        particle.acceleration = lorentz

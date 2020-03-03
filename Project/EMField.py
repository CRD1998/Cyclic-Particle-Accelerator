import log
import numpy as np
import scipy.constants as const
import math

class EMField:
    """
    Class to represent an  electromagnetic field.

    Parse in the electric field at its maximum and it will be converted to a time-varying
    electric field of the form Acos(wt+phi) 

    Electric field in N/C or V/m
    Magnetic field in T
    """

    def __init__(self, ElectricField = np.array([0,0,0], dtype=float), 
                MagneticField = np.array([0,0,0], dtype=float), 
                ElectricFieldWidth = np.array([-1,1], dtype=float)):
        self.magnetic = np.array(MagneticField,dtype=float)
        self.electric = np.array(ElectricField,dtype=float)
        self.electricLowerBound  = float(ElectricFieldWidth[0])
        self.electricUpperBound  = float(ElectricFieldWidth[1])
        log.logger.info('electromagentic field generated')

    def __repr__(self):
        return 'EM Field: \n Electric Field = {0}  \n Magnetic Field = {1}'.format(self.electric,self.magnetic)

    def electricMag(self):
        return np.linalg.norm(self.electric)
    
    def magneticMag(self):
        return np.linalg.norm(self.magnetic)
    
    def frequency(self,particle):
        return abs(particle.charge)*self.magneticMag()/particle.mass

    def getAcceleration(self,particleBunch,time,deltaT):
        for particle in particleBunch:
            electricField = lambda x: [i*math.cos(self.frequency(particle)*time) for i in self.electric] if self.electricLowerBound < x < self.electricUpperBound else [0,0,0]
            lorentz = np.array(electricField(particle.position[0]),dtype=float)
            lorentz += np.cross(particle.velocity, self.magnetic)
            lorentz *= particle.charge/ (particle.mass*particle.gamma())
            particle.acceleration = lorentz
            particle.eulerCromer(deltaT)

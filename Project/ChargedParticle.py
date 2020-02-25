import numpy as np
from Particle import Particle
import scipy.constants as const
import logging

class ChargedParticle(Particle):
    """
    Class to simulate a charged particle, it is a child of the Particle class in Particle.py

    Charge in Coulombs
    """
    
    def __init__(self,  Name='Point', Mass=1.0, Charge=1.0, Position=np.array([0,0,0], dtype=float), Velocity=np.array([0,0,0], dtype=float), Acceleration=np.array([0,0,0], dtype=float)):
        super().__init__(Name=Name, Mass=Mass, Position=Position, Velocity=Velocity, Acceleration=Acceleration)
        self.charge = float(Charge)
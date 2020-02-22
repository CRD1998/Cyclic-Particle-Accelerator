import numpy as np 
import matplotlib.pyplot as plt 
import random
from ChargedParticle import ChargedParticle

"""
This script will not form part of the final simulation, this script is just providing the code that
will be used in the bunch class. From this, it seems that the bunch class can just take the number
of particles in the bunch and the desired average position, velcoity and acceleration. Methods,
similar to the functions below can then be used to generate random vectors. 

It may also take the desired intial spread.

#TODO look into avoiding repeated vectors
"""

n = 10
average_position = np.array([10,20,30])

def float_random_coordinates(n, average_position):
    x_positions = np.array([random.uniform(average_position[0]*0.90, average_position[0]*1.10) for _ in range(n)])
    y_positions = np.array([random.uniform(average_position[1]*0.90, average_position[1]*1.10) for _ in range(n)])
    z_positions = np.array([random.uniform(average_position[2]*0.90, average_position[2]*1.10) for _ in range(n)])
    return x_positions,y_positions,z_positions

def integer_random_coordinates(n, average_position):
    x_positions = np.array([random.randint(average_position[0]*0.90, average_position[0]*1.10) for _ in range(n)])
    y_positions = np.array([random.randint(average_position[1]*0.90, average_position[1]*1.10) for _ in range(n)])
    z_positions = np.array([random.randint(average_position[2]*0.90, average_position[2]*1.10) for _ in range(n)])
    return x_positions,y_positions,z_positions

x, y, z = integer_random_coordinates(n, average_position)
particle_bunch, no = [], 1
for (i,j,k) in zip(x,y,z):
    proton = ChargedParticle('proton-'+str(no),1,1,[i,j,k])
    no += 1
    particle_bunch.append(proton)
for i in particle_bunch:
    print(i.name, i.position)
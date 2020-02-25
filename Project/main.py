import log
import numpy as np
from Particle import Particle
from ChargedParticle import ChargedParticle
from EMField import EMField
import scipy.constants as const
import matplotlib.pyplot as plt
import copy

proton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [0,3000,0])
field = EMField([0,0.1,0], [0,0,1.6*10**(-5)])

time, deltaT, duration = 0, 10**(-5), 0.03

Data = []

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    field.getAcceleration(proton)
    proton.eulerCromer(deltaT)
    temp_particle= copy.deepcopy(proton)
    Data.append(temp_particle.position)
log.logger.info('simulation finished')

log.logger.info('creating plots')
x, y = [i[0] for i in Data], [i[1] for i in Data]
magneticX, magneticY = np.meshgrid(list(range(-1,6)), list(range(-3,4)))

plt.figure('Proton in EM Field')
plt.plot(x,y,color='blue',label="Euler Cromer Path")

plt.xlabel(r'x  position  [$m$]')
plt.ylabel(r'y  position [$m$]')
plt.title('Plot to show the path of a proton in an electromagnetic field')
plt.legend()
plt.show()
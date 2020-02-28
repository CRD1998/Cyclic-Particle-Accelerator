import log
import numpy as np
from ChargedParticle import ChargedParticle
from Bunch import Bunch
from EMField import EMField
import scipy.constants as const
import matplotlib.pyplot as plt
import copy

field = EMField([.1,0,0], [0,0,1.6*10**(-5)])
proton = ChargedParticle('proton',const.m_p,const.e)
protons = Bunch(proton,0.047,3)

time, deltaT, duration = 0, 10**(-5), 0.03

Data = []
for i in protons.bunch:
    Data.append([])

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    for (particle,array) in zip(protons.bunch,Data):
        field.getAcceleration(particle,time)
        particle.eulerCromer(deltaT)
        temp_particle = copy.deepcopy(particle)
        array.append(temp_particle)
log.logger.info('simulation finished')

magneticX, magneticY = np.meshgrid(list(range(-20,21,5)), list(range(-20,16,5)))

log.logger.info('creating plots')
plt.figure('Cyclotron with Bunch')
for i in Data:
    x,y = [j.position[0] for j in i], [j.position[1] for j in i]
    final = [x[-1], y[-1]]
    plt.plot(x,y,label=i[0].name)
    plt.scatter(final[0], final[1], marker='.')

plt.axvspan(field.electricLowerBound,field.electricUpperBound,alpha=0.5, color='grey',label='Electric Field')
plt.scatter(magneticX,magneticY,marker=r'$\bigotimes$',s=95,color='black',label='Magnetic Field')
plt.xlabel(r'x  position  [m]')
plt.ylabel(r'y  position [m]')
plt.legend(loc='upper right',framealpha=1)
plt.show()
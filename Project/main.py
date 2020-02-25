import log
import numpy as np
from Particle import Particle
from ChargedParticle import ChargedParticle
from EMField import EMField
import scipy.constants as const
import matplotlib.pyplot as plt
import copy

proton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [3000,0,0])
field = EMField([0.1,0,0], [0,0,1.6*10**(-5)])

time, deltaT, duration = 0, 10**(-5), 0.03

Data = []
electricField = []
timeSeries = []

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    timeSeries.append(time)
    field.getAcceleration(proton,time)
    proton.eulerCromer(deltaT)
    temp_particle, temp_field= copy.deepcopy(proton), copy.deepcopy(field.ImplementElectricField(proton,time))
    Data.append(temp_particle.position)
    electricField.append(temp_field)
log.logger.info('simulation finished')

log.logger.info('creating plots')
x, y = [i[0] for i in Data], [i[1] for i in Data]
final = [x[-1], y[-1]]
magneticX, magneticY = np.meshgrid(list(range(-15,21,5)), list(range(-20,16,5)))

#caption = r' For this cyclotron: $\vec{B}=(1.6\times 10^{-5})\mathbf{\hat{e}_z}$ $T$ , $\vec{E}=0.1\mathbf{\hat{e}_y}$ $NC^{-1}$'

plt.figure('Proton in EM Field')
#plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment='center', fontsize=10)
plt.plot(x,y,color='blue',label="Euler Cromer Path")
plt.scatter(final[0],final[1],color='blue',marker='o')
plt.xlabel(r'x  position  [$m$]')
plt.ylabel(r'y  position [$m$]')
plt.title('Plot to show the path of a proton in an electromagnetic field')
plt.axvspan(-0.98,0.98,alpha=0.5, color='grey',label='Electric Field')
plt.scatter(magneticX,magneticY,marker=r'$\bigotimes$',s=95,color='black',label='Magnetic Field')
plt.legend(loc='upper right',framealpha=1)

plt.figure('Electric Field')
plt.plot(timeSeries,[i[0] for i in electricField],label='x component of electric field')
plt.xlabel('time [$s$]')
plt.ylabel('Electric Field [$NC^{-1}$]')
plt.title(r'Plot to show $E_{x}$ of $\vec{E}$')

plt.show()
import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([.1,0,0], [0,0,1.6*10**(-5)])
protons = ProtonBunch(0.047)

log.logger.info('Initial average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Initial average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Initial bunch spread: %s m' % protons.spread())

time, deltaT, duration = 0, 10**(-6), 0.005

timeSeries = []
Data = []

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    timeSeries.append(time)
    field.getAcceleration(protons.bunch, time, deltaT)
    temp_bunch = copy.deepcopy(protons)
    Data.append(temp_bunch)
log.logger.info('simulation finished')

#magneticX, magneticY = np.meshgrid(list(range(-20,21,5)), list(range(-20,16,5)))

log.logger.info('creating plots')
plt.figure('Cyclotron with Bunch')

plt.axvspan(field.electricLowerBound,field.electricUpperBound,alpha=0.5, color='grey',label='Electric Field')
#plt.scatter(magneticX,magneticY,marker=r'$\bigotimes$',s=95,color='black',label='Magnetic Field')

x,y,y_downSpread,y_upSpread = [],[],[],[]
for bunch in Data:
    x.append(bunch.averagePosition()[0])
    y.append(bunch.averagePosition()[1])
    y_downSpread.append(bunch.averagePosition()[1] - abs(bunch.spread()[1]))
    y_upSpread.append(bunch.averagePosition()[1] + abs(bunch.spread()[1]))
final = [x[-1], y[-1]]
plt.plot(x,y,label=protons.bunchName)
plt.scatter(final[0], final[1], marker='.')
plt.fill_between(x,y_downSpread,y_upSpread)
plt.xlabel(r'x  position  [m]')
plt.ylabel(r'y  position [m]')
plt.legend(loc='upper left',framealpha=1)



plt.show()
import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

#field = EMField([10*8,0,0], [0,0,0.006568], [-0.25,0.25], 7*np.pi/4) 
#protons = ProtonBunch(51511,2)
field = EMField([50*10**3,0,0], [0,0,2.82], [-0.05,0.05], 7*np.pi/4) 
protons = ProtonBunch(120*10**(6),1)

log.logger.info('Initial average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Initial average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Initial average position %s m' % protons.averagePosition())
log.logger.info('Initial bunch position spread: %s m' % protons.positionSpread())
log.logger.info('Initial bunch energy spread: %s eV' % protons.energySpread())

method = 2
time, deltaT, duration = 0, 2*10**(-11), 2.38*10**(-8)*100

inital_bunch = copy.deepcopy(protons)

timeSeries = [0.]
Data = [inital_bunch]

log.logger.info('starting simulation')
while time <= duration:
    dt = protons.adaptiveStep(deltaT,field)
    time += dt
    timeSeries.append(time)
    field.getAcceleration(protons.bunch, time, dt)
    protons.update(dt,field,time,method)
    temp_bunch = copy.deepcopy(protons)
    Data.append(temp_bunch)

log.logger.info('simulation finished')

log.logger.info('Final average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Final average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Final average position %s m' % protons.averagePosition())
log.logger.info('Final bunch position spread: %s m' % protons.positionSpread())
log.logger.info('Final bunch energy spread: %s eV' % protons.energySpread())

log.logger.info('writing lists to file')
np.savez('synchrocyclotron_data', time=timeSeries, protons_data=Data)
log.logger.info('file writing complete')

""""
log.logger.info('building lists')

x,y = [],[]
for bunch in Data:
    x.append(bunch.averagePosition()[0])
    y.append(bunch.averagePosition()[1])

final = [x[-1], y[-1]]
magneticX, magneticY = np.meshgrid(list(range(-6,7,2)), list(range(-10,1,2)))

log.logger.info('creating plots')

plt.figure('Synchrocyclotron Bunch with Spread')
plt.axvspan(field.electricLowerBound,field.electricUpperBound,alpha=0.5, color='grey',label='Electric Field',zorder=1)
#plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.plot(x,y,label=protons.bunchName,color='blue',zorder=3)
plt.scatter(final[0], final[1],color='blue',zorder=3)
plt.xlabel(r'x  position  [m]')
plt.ylabel(r'y  position [m]')
plt.legend(loc='upper left',framealpha=1)

plt.show()
"""
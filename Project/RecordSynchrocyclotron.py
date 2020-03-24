import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([10*8,0,0], [0,0,0.006568], [-0.25,0.25], 7*np.pi/4) 
protons = ProtonBunch(51511,2)
#protons_2 = ProtonBunch(0.0047,25)
#protons = protons_1 + protons_2

log.logger.info('Initial average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Initial average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Initial average position %s m' % protons.averagePosition())
log.logger.info('Initial bunch position spread: %s m' % protons.positionSpread())
log.logger.info('Initial bunch energy spread: %s eV' % protons.energySpread())

time, deltaT, duration = 0, 10**(-8), 10**(-5)*5

inital_bunch = copy.deepcopy(protons)

timeSeries = [0.]
Data = [inital_bunch]

log.logger.info('starting simulation')
while time <= duration:
    dt = protons.adaptiveStep(deltaT,field)
    time += dt
    timeSeries.append(time)
    field.getAcceleration(protons.bunch, time, dt)
    protons.update(dt,field,time,2)
    temp_bunch = copy.deepcopy(protons)
    Data.append(temp_bunch)

log.logger.info('simulation finished')

log.logger.info('Final average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Final average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Final average position %s m' % protons.averagePosition())
log.logger.info('Final bunch position spread: %s m' % protons.positionSpread())
log.logger.info('Final bunch energy spread: %s eV' % protons.energySpread())

log.logger.info('building lists')

x,y,z,x_upSpread,x_downSpread = [],[],[],[],[]
for bunch in Data:
    x.append(bunch.averagePosition()[0])
    y.append(bunch.averagePosition()[1])
    z.append(bunch.averagePosition()[2])

    x_downSpread.append(bunch.averagePosition()[0] - bunch.positionSpread()[0])
    x_upSpread.append(bunch.averagePosition()[0] + bunch.positionSpread()[0])

final = [x[-1], y[-1]]
magneticX, magneticY = np.meshgrid(list(range(-3,4,1)), list(range(-3,3,1)))

log.logger.info('creating plots')

plt.figure('Cyclotron Bunch with Spread')
plt.axvspan(field.electricLowerBound,field.electricUpperBound,alpha=0.5, color='grey',label='Electric Field',zorder=1)
#plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field')
plt.plot(x,y,label=protons.bunchName,color='blue')
plt.fill_betweenx(y, x_downSpread, x_upSpread, color='cornflowerblue')
plt.scatter(final[0], final[1],color='blue')
plt.xlabel(r'x  position  [m]')
plt.ylabel(r'y  position [m]')
plt.legend(loc='upper left',framealpha=1)

plt.show()
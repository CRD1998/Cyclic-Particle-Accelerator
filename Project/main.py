import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([.1,0,0], [0,0,1.6*10**(-5)], [-0.026,0.026]) 
protons = ProtonBunch(0.047)

log.logger.info('Initial average kinetic energy: %s eV' % protons.KineticEnergy())
log.logger.info('Initial average momentum: %s kg m/s' % protons.momentum())
log.logger.info('Initial average position %s m' % protons.averagePosition())
log.logger.info('Initial bunch spread: %s m' % protons.spread())

time, deltaT, duration = 0, 10**(-5), 0.0042*5

inital_bunch = copy.deepcopy(protons)

timeSeries = [0.]
Data = [inital_bunch]

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    timeSeries.append(time)
    field.getAcceleration(protons.bunch, time, deltaT)
    temp_bunch = copy.deepcopy(protons)
    Data.append(temp_bunch)

log.logger.info('simulation finished')
log.logger.info('building lists')

x,y,y_downSpread,x_upSpread,x_downSpread,y_upSpread = [],[],[],[],[],[]
bunch_spread_x, bunch_spread_y = [], []
for bunch in Data:
    x.append(bunch.averagePosition()[0])
    y.append(bunch.averagePosition()[1])

    x_downSpread.append(bunch.averagePosition()[0] - bunch.spread()[0])
    x_upSpread.append(bunch.averagePosition()[0] + bunch.spread()[0])
    y_downSpread.append(bunch.averagePosition()[1] - bunch.spread()[1])
    y_upSpread.append(bunch.averagePosition()[1] + bunch.spread()[1])

    bunch_spread_x.append(bunch.spread()[0])
    bunch_spread_y.append(bunch.spread()[1])

final = [x[-1], y[-1]]
magneticX, magneticY = np.meshgrid(list(range(-4,4,1)), list(range(-5,5,1)))


log.logger.info('creating plots')

plt.figure('Cyclotron Bunch with Spread')
plt.axvspan(field.electricLowerBound,field.electricUpperBound,alpha=0.5, color='grey',label='Electric Field')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field')
plt.plot(x,y,label=protons.bunchName,color='blue')
plt.fill_between(x,y_downSpread,y_upSpread, color='cornflowerblue',label="spread")
plt.fill_betweenx(y, x_downSpread, x_upSpread, color='cornflowerblue')
plt.scatter(final[0], final[1],color='blue')
plt.xlabel(r'x  position  [m]')
plt.ylabel(r'y  position [m]')
plt.legend(loc='upper left',framealpha=1)

plt.figure('Spread')
plt.plot(timeSeries, bunch_spread_x, label=r'$\sigma_{x}$')
plt.plot(timeSeries, bunch_spread_y, label=r'$\sigma_{y}$')
plt.xlabel('time [s]')
plt.ylabel('spread [m]')
plt.legend(loc='best',framealpha=1)

plt.show()
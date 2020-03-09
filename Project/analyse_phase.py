import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([.1,0,0], [0,0,1.6*10**(-5)], [-0.026,0.026]) 
protons = ProtonBunch(0.0047,20)
inital_bunch = copy.deepcopy(protons)

phases = [0, np.pi/4, np.pi/2]

deltaT, duration = 10**(-5), 0.0042*20
timeSeries = []
Data = []
for i in range(len(phases)):
    Data.append([])

log.logger.info('starting simulation')

for (angle,loop) in zip(phases,Data):
    time = 0
    field.phase = angle
    protons = copy.deepcopy(inital_bunch)
    log.logger.info('phase currently being investigated: %s radians' % angle)
    while time <= duration:
        time += deltaT
        if angle == phases[0]: # only save data to the time series on the final iteration
            timeSeries.append(time)
        field.getAcceleration(protons.bunch, time, deltaT)
        temp_spread = copy.deepcopy(protons.spread())
        loop.append(temp_spread)

log.logger.info('simulation finished')

plt.figure('Phase Spreads')
for (i,phase) in zip(Data,phases):
    x_spread = [j[1] for j in i]
    plt.plot(timeSeries,x_spread,label='Phase: '+str(round(phase,2)))
plt.xlabel('time [s]')
plt.ylabel(r'$\sigma_x$ [m]')
plt.legend()

plt.figure('Linear Spreads')
for (i,phase) in zip(Data,phases):
    x_spread = [j[0] for j in i]
    linear_fit = np.polynomial.polynomial.Polynomial.fit(timeSeries,x_spread,1)
    line = np.poly1d(linear_fit.coef)
    spread = [line(t) for t in timeSeries]
    plt.plot(timeSeries,spread,label='Phase: '+str(round(phase,2)))
plt.xlabel('time [s]')
plt.ylabel(r'$\sigma_x$ [m]')
plt.legend()

plt.show()
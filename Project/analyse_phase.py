import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes,InsetPosition,mark_inset)
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([.1,0,0], [0,0,1.6*10**(-5)], [-0.026,0.026]) 
protons = ProtonBunch(0.0047,100)
inital_bunch = copy.deepcopy(protons)

phases = [i*np.pi/8 for i in range(0,16)]
phases = [angle for angle in phases if np.cos(angle) >= -0.00001]

deltaT, duration = 10**(-5), 0.0042*5
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

fig, ax = plt.subplots()
ax1 = plt.axes([0,0,1,1]) # create inset instance
ip = InsetPosition(ax1, [0.2,0.135,0.55,0.55]) # all of these are fractional, format: x coordinate of plot. y coordinate of plot, height, width
ax1.set_axes_locator(ip) # assign inset position to inset instance
mark_inset(ax, ax1, loc1=3, loc2=4, fc="none", ec='0.5') # draw grey lines from inset position to data on plot
for (i,phase) in zip(Data,phases):
    x_spread = [j[0] for j in i]
    linear_fit = np.polynomial.polynomial.Polynomial.fit(timeSeries,x_spread,1)
    line = np.poly1d(linear_fit.coef)
    spread = [line(t) for t in timeSeries]
    ax.plot(timeSeries,spread,label='Phase: '+str(round(phase,2)))
    ax1_x, ax1_y = [], []
    for (x,y) in zip(timeSeries,spread):
        if 0.<y<0.00005: # the data we are focusing on
            ax1_x.append(x)
            ax1_y.append(y)
    ax1.plot(ax1_x,ax1_y) # plot the data on the inset
ax.set_xlabel('time [s]')
ax.set_ylabel(r'$\sigma_x$ [m]')
ax.legend()
plt.savefig('Linear_Spreads_Zoom.png')

plt.figure('Phases')
theta = np.linspace(0,2*np.pi,10000)
cosine = [np.cos(x) for x in theta]
phase_points = [np.cos(x) for x in phases]
plt.plot(theta, cosine, label=r'$cos(\theta)$', color='red',)
plt.scatter(phases, phase_points, label='Electric field phases',color='black')
plt.xlabel(r'$\theta$ (radians)')
plt.legend()

plt.show()
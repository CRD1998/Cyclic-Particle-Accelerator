import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

"""
This file will analyse several different electric field phase shifts. It will then plot the spread of
the x-position of a proton bunch as well as the bunch's spread of kinetic energy against time. It will
then fit a linear curve to the position and energy spreads so the variation of the spreads against time
can be more easily visualised.

The 'phase_data.npz' contains the position and kinetic energy spreads for a proton bunch containing
100 protons over approximately 100 revolutions. If you do not have the file in the same location as
this file, the simulation will run and generate it for you. I would recommend you reduce the size of the
bunch and the simulation's duration as running the simulation for 100 revolutions with a 100 particle
bunch size took just over six hours to run on my machine.
"""

phase_keys = [
    r'$0$', r'$\frac{1}{8} \pi$', r'$\frac{1}{4} \pi$', r'$\frac{3}{8} \pi$',
    r'$\frac{1}{2} \pi$', r'$\frac{5}{8} \pi$', r'$\frac{3}{4} \pi$', r'$\frac{7}{8} \pi$',
    r'$\pi$', r'$\frac{9}{8} \pi$', r'$\frac{5}{4} \pi$', r'$\frac{11}{8} \pi$',
    r'$\frac{3}{2} \pi$', r'$\frac{13}{8} \pi$', r'$\frac{7}{4} \pi$', r'$\frac{15}{8} \pi$'
]
phase_values = [i*np.pi/8 for i in range(0,16)]
phase_dict = {phase_keys[i]:phase_values[i] for i in range(len(phase_keys)) if np.cos(phase_values[i])>=-0.00001}
phases = list(phase_dict.values())

def generate_file(phases):
    field = EMField([.1,0,0], [0,0,1.6*10**(-5)], [-0.026,0.026]) 
    protons = ProtonBunch(0.0047,100)
    inital_bunch = copy.deepcopy(protons)

    deltaT, duration = 10**(-5), 0.0041*101
    timeSeries = []
    positionSpread = []
    energySpread = []
    for _ in range(len(phases)):
        positionSpread.append([])
        energySpread.append([])

    log.logger.info('starting simulation')

    for (angle,loop_1,loop_2) in zip(phases,positionSpread,energySpread):
        time = 0
        field.phase = angle
        protons = copy.deepcopy(inital_bunch)
        log.logger.info('phase currently being investigated: %s radians' % angle)
        while time <= duration:
            time += deltaT
            if angle == phases[0]: # only save data to the time series on the first iteration
                timeSeries.append(time)
            field.getAcceleration(protons.bunch, time, deltaT)
            protons.update(deltaT,field,time)
            temp_spread = copy.deepcopy(protons.positionSpread())
            temp_energy = copy.deepcopy(protons.energySpread())
            loop_1.append(temp_spread)
            loop_2.append(temp_energy/const.e)
    log.logger.info('simulation finished')

    log.logger.info('writing to file')
    np.savez('phase_data', time=timeSeries, positions=positionSpread, energies=energySpread)
    log.logger.info('file written')
    return np.load('phase_data.npz')

try:
    simulation_data = np.load('phase_data.npz')
except FileNotFoundError:
    simulation_data = generate_file(phases)

timeSeries = simulation_data['time']
positions = simulation_data['positions']
energies = simulation_data['energies']

plt.figure('Phase Position Spreads')
for (i,phase) in zip(positions,phase_dict.keys()):
    x_spread = [j[1] for j in i]
    plt.plot(timeSeries,x_spread,label='Phase: '+ phase)
plt.xlabel('time [s]')
plt.ylabel(r'$\sigma_x$ [m]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Phase Energy Spreads')
for (i,phase) in zip(energies,phase_dict.keys()):
    plt.plot(timeSeries,i,label='Phase: '+ phase)
plt.xlabel('time [s]')
plt.ylabel(r'Kinetic Energy Spread ($1\sigma$) [eV]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Linear Position Phase Spreads')
for (i,phase) in zip(positions,phase_dict.keys()):
    x_spread = [j[0] for j in i]
    linear_fit = np.polynomial.polynomial.Polynomial.fit(timeSeries,x_spread,1)
    line = np.poly1d(linear_fit.coef)
    spread = [line(t) for t in timeSeries]
    plt.plot(timeSeries,spread,label='Phase: '+phase)
plt.xlabel('time [s]')
plt.ylabel(r'$\sigma_x$ [m]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Linear Energy Phase Spreads')
for (i,phase) in zip(energies,phase_dict.keys()):
    linear_fit = np.polynomial.polynomial.Polynomial.fit(timeSeries,i,1)
    line = np.poly1d(linear_fit.coef)
    spread = [line(t) for t in timeSeries]
    plt.plot(timeSeries,spread,label='Phase: '+phase)
plt.xlabel('time [s]')
plt.ylabel(r'Kinetic Energy Spread ($1\sigma$) [eV]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Phases')
theta = np.linspace(0,2*np.pi,10000)
cosine = [np.cos(x) for x in theta]
phase_points = [np.cos(x) for x in phases]
plt.plot(theta, cosine, label=r'$cos(\theta)$', color='red',zorder=1)
plt.scatter(phases, phase_points, label='Electric field phases',color='black',zorder=10)
for x,y,name in zip(phases,phase_points,phase_dict.keys()):
    if phases.index(x) < 5:
        plt.text(x+0.15,y+0.1,name)
    else:
        plt.text(x-0.25,y+0.1,name)
plt.ylim(-1.25,1.25)
plt.xlabel(r'$\theta$ (radians)')
plt.legend(loc='lower left')
plt.tick_params(which='both',direction='in',right=True,top=True)

# format plots so that they look similar to QtiPlot
matplotlib.rcParams['lines.linewidth'] = 6
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.major.size'] = 9
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 1.9
matplotlib.rcParams['xtick.minor.width'] = 1.3
matplotlib.rcParams['ytick.major.size'] = 9
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.width'] = 1.9
matplotlib.rcParams['ytick.minor.width'] = 1.3

plt.show()
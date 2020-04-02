import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes,InsetPosition,mark_inset)
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

"""
foo
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
    field = EMField([500000,0,0], [0,0,2.82], [-0.05,0.05]) 
    protons = ProtonBunch(120*10**(6),25)
    inital_bunch = copy.deepcopy(protons)

    deltaT, duration = 2*10**(-10), 2.78*10**(-8)*10
    timeSeries = []
    positionSpread = []
    energySpread = []
    for _ in range(len(phases)):
        positionSpread.append([])
        energySpread.append([])
        timeSeries.append([])

    log.logger.info('starting simulation')

    for (angle,loop_1,loop_2,loop_3) in zip(phases,positionSpread,energySpread,timeSeries):
        time = 0
        field.phase = angle
        protons = copy.deepcopy(inital_bunch)
        log.logger.info('phase currently being investigated: %s radians' % angle)
        while time <= duration:
            dt = protons.adaptiveStep(deltaT,field)
            time += dt
            field.getAcceleration(protons.bunch, time, dt)
            protons.update(dt,field,time,2)
            temp_spread = copy.deepcopy(protons.positionSpread())
            temp_energy = copy.deepcopy(protons.energySpread())
            loop_1.append(temp_spread)
            loop_2.append(temp_energy)
            loop_3.append(time)
    log.logger.info('simulation finished')

    log.logger.info('writing to file')
    np.savez('phase_synchro_data', time=timeSeries, positions=positionSpread, energies=energySpread)
    log.logger.info('file written')
    return np.load('phase_synchro_data.npz',allow_pickle=True)

try:
    simulation_data = np.load('phase_synchro_data.npz',allow_pickle=True)
except FileNotFoundError:
    simulation_data = generate_file(phases)

times = simulation_data['time']
positions = simulation_data['positions']
energies = simulation_data['energies']
timeSeries = []
for series in times:
    timeSeries.append([])
    for i in series:
        timeSeries[len(timeSeries)-1].append(i*10**(6))

plt.figure('Synchro Phase Position Spreads')
for (i,phase,time) in zip(positions,phase_dict.keys(),timeSeries):
    y_spread = [j[1] for j in i]
    plt.plot(time,y_spread,label='Phase: '+ phase)
plt.xlabel(r'time [$\mu$s]')
plt.ylabel(r'$\sigma_y$ [m]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Synchro Phase Energy Spreads')
for (i,phase,time) in zip(energies,phase_dict.keys(),timeSeries):
    plt.plot(time,[j*10**(-6) for j in i],label='Phase: '+ phase)
plt.xlabel(r'time [$\mu$s]')
plt.ylabel(r'Kinetic Energy Spread ($1\sigma$) [MeV]')
plt.legend(loc='upper left')
plt.tick_params(which='both',direction='in',right=True,top=True)

fig, ax = plt.subplots()
for (i,phase,time) in zip(energies,phase_dict.keys(),timeSeries):
    adjusted_energies = [j*10**(-6) for j in i]
    linear_fit = np.polynomial.polynomial.Polynomial.fit(time,adjusted_energies,1)
    line = np.poly1d(linear_fit.coef)
    a, b = line[0], line[1]
    spread = [a*t+b for t in time]
    ax.plot(time,spread,label='Phase: '+phase)
ax.set_xlabel(r'time [$\mu$s]')
ax.set_ylabel(r'Kinetic Energy Spread ($1\sigma$) [MeV]')
ax.legend(loc='upper left')
ax.tick_params(which='both',direction='in',right=True,top=True)

plt.show()
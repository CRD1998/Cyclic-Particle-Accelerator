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
timesteps = [10**(-6), 10**(-5), 10**(-4)]
step_strings = ['1e-6', '1e-5', '1e-4']
step_dict = {step_strings[i]:timesteps[i] for i in range(len(timesteps))}

def generate_file():
    field = EMField([0,0,0], [0,0,1.6*10**(-5)], [0,0]) 
    protons = ProtonBunch(0.0047,100)
    inital_bunch = copy.deepcopy(protons)

    duration = 0.0041*100
    timeSeries = []
    kineticEnergies = []
    momenta = []
    timeSeries = []
    for _ in range(len(timesteps)):
        kineticEnergies.append([])
        momenta.append([])
        timeSeries.append([])

    log.logger.info('starting simulation')

    for (step,loop_1,loop_2,loop_3) in zip(timesteps,timeSeries,kineticEnergies,momenta):
        time = 0
        deltaT = step
        field.phase = step
        protons = copy.deepcopy(inital_bunch)
        log.logger.info('phase currently being investigated: %s radians' % step)
        while time <= duration:
            time += deltaT
            field.getAcceleration(protons.bunch, time, deltaT)
            protons.update(deltaT,field,time,2)
            temp_energy = copy.deepcopy(protons.KineticEnergy())
            temp_momentum = copy.deepcopy(protons.momentum())
            loop_1.append(time)
            loop_2.append(temp_energy)
            loop_3.append(np.linalg.norm(temp_momentum))
    log.logger.info('simulation finished')

    log.logger.info('writing to file')
    np.savez('timestep_data', time=timeSeries, energies=kineticEnergies, momenta=momenta)
    log.logger.info('file written')
    return np.load('timestep_data.npz',allow_pickle=True)

try:
    simulation_data = np.load('timestep_data.npz',allow_pickle=True)
except FileNotFoundError:
    simulation_data = generate_file()

timeSeries = simulation_data['time']
energies = simulation_data['energies']
momenta = simulation_data['momenta']

fig, ax = plt.subplots()
ax1 = plt.axes([0,0,1,1]) # create inset instance
ip = InsetPosition(ax1, [0.4,0.4,0.45,0.45]) # all of these are fractional, format: x coordinate of plot. y coordinate of plot, height, width
ax1.set_axes_locator(ip) # assign inset position to inset instance
mark_inset(ax, ax1, loc1=1, loc2=2, fc="none", ec='0.5') # draw grey lines from inset position to data on plot
for i,time,step in zip(momenta,timeSeries,step_dict.keys()):
    fractional_momenta = [j/i[0] for j in i]
    ax.plot(time, fractional_momenta,label=step)
    if step == '1e-4':
        continue
    ax1_x, ax1_y = [],[]
    for x,y in zip(time,fractional_momenta):
        if 0.<x<0.4: # the data we are focusing on
            ax1_x.append(x)
            ax1_y.append(y)
    ax1.plot(ax1_x,ax1_y) # plot the data on the inset
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'Fractional $\|\vec{p}\|$')
ax.legend(loc='lower left')
ax.tick_params(which='both',direction='in',right=True,top=True)
plt.savefig('timestep-momentum-constantB.png')

fig, ax = plt.subplots()
ax1 = plt.axes([0,0,1,1]) # create inset instance
ip = InsetPosition(ax1, [0.4,0.4,0.45,0.45]) # all of these are fractional, format: x coordinate of plot. y coordinate of plot, height, width
ax1.set_axes_locator(ip) # assign inset position to inset instance
mark_inset(ax, ax1, loc1=1, loc2=2, fc="none", ec='0.5') # draw grey lines from inset position to data on plot
for i,time,step in zip(energies,timeSeries,step_dict.keys()):
    fractional_energies = [j/i[0] for j in i]
    ax.plot(time, fractional_energies,label=step)
    if step == '1e-4':
        continue
    ax1_x, ax1_y = [],[]
    for x,y in zip(time,fractional_energies):
        if 0.<x<0.4: # the data we are focusing on
            ax1_x.append(x)
            ax1_y.append(y)
    ax1.plot(ax1_x,ax1_y) # plot the data on the inset
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'Fractional $E_k$')
ax.legend(loc='lower left')
ax.tick_params(which='both',direction='in',right=True,top=True)
plt.savefig('timestep-energies-constantB.png')
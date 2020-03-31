import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

try:
    simulation_data = np.load('methods_data.npz', allow_pickle=True)
except FileNotFoundError:
    try:
        simulation_data = np.load('cyclotron_data.npz', allow_pickle=True)
    except FileNotFoundError:
        import RecordCyclotron
        simulation_data = np.load('cyclotron_data.npz', allow_pickle=True)

time = simulation_data['time']
eulerData = simulation_data['euler']
cromerData = simulation_data['cromer']
verletData = simulation_data['verlet']
RK4Data = simulation_data['rk']

fractional_kinetic_euler = [i.KineticEnergy()/eulerData[0].KineticEnergy() for i in eulerData]
fractional_kinetic_cromer = [i.KineticEnergy()/cromerData[0].KineticEnergy() for i in cromerData]
fractional_kinetic_verlet = [i.KineticEnergy()/verletData[0].KineticEnergy() for i in verletData]
fractional_kinetic_rk4 = [i.KineticEnergy()/RK4Data[0].KineticEnergy() for i in RK4Data]

fractional_momentum_euler = [np.linalg.norm(i.momentum())/np.linalg.norm(eulerData[0].momentum()) for i in eulerData]
fractional_momentum_cromer = [np.linalg.norm(i.momentum())/np.linalg.norm(cromerData[0].momentum()) for i in cromerData]
fractional_momentum_verlet = [np.linalg.norm(i.momentum())/np.linalg.norm(verletData[0].momentum()) for i in verletData]
fractional_momentum_rk4 = [np.linalg.norm(i.momentum())/np.linalg.norm(RK4Data[0].momentum()) for i in RK4Data]

plt.figure('Euler-EulerCromer Fractional Kinetic Energy')
plt.plot(time, fractional_kinetic_euler, label='Euler')
plt.plot(time, fractional_kinetic_cromer, label='Euler-Cromer')
plt.xlabel('time [s]')
plt.ylabel('Fractional Kinetic Energy')
plt.legend()
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Euler-EulerCromer Fractional Momentum')
plt.plot(time, fractional_momentum_euler, label='Euler')
plt.plot(time, fractional_momentum_cromer, label='Euler-Cromer')
plt.xlabel('time [s]')
plt.ylabel(r'Fractional $\parallel\vec{p}\parallel$')
plt.legend()
plt.tick_params(which='both',direction='in',right=True,top=True)


plt.figure('Verlet-RK4 Fractional Kinetic Energy')
plt.plot(time, fractional_kinetic_verlet, label='Velocity Verlet')
plt.plot(time, fractional_kinetic_rk4, label='4th Order RK')
plt.xlabel('time [s]')
plt.ylabel('Fractional Kinetic Energy')
plt.legend()
plt.tick_params(which='both',direction='in',right=True,top=True)
plt.ticklabel_format(useOffset=False)

plt.figure('Verlet-RK4 Fractional Momentum')
plt.plot(time, fractional_momentum_verlet, label='Velocity Verlet')
plt.plot(time, fractional_momentum_rk4, label='4th Order RK')
plt.xlabel('time [s]')
plt.ylabel(r'Fractional $\parallel\vec{p}\parallel$')
plt.legend()
plt.tick_params(which='both',direction='in',right=True,top=True)
plt.ticklabel_format(useOffset=False)

plt.show()
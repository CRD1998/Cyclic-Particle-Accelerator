import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib

try:
    simulation_data = np.load('cyclotron_data.npz',allow_pickle=True)
except FileNotFoundError:
    import RecordCyclotron
    simulation_data = np.load('cyclotron_data.npz',allow_pickle=True)

try:
    sychrocyclotron = np.load('synchrocyclotron_data.npz', allow_pickle=True)
except FileNotFoundError:
    import RecordSynchrocyclotron
    sychrocyclotron = np.load('synchrocyclotron_data.npz', allow_pickle=True)

eulerData = simulation_data['euler']
cromerData = simulation_data['cromer']
verletData = simulation_data['verlet']
RK4Data = simulation_data['rk']

data = sychrocyclotron['protons_data']

euler_x, euler_y = [], []
cromer_x, cromer_y = [], []
verlet_x, verlet_y = [], []
rk4_x, rk4_y = [], []
synchro_x , synchro_y = [], []
for i,j,k,l in zip(eulerData,cromerData,verletData,RK4Data):
    position = i.averagePosition()
    euler_x.append(position[0])
    euler_y.append(position[1])

    position = j.averagePosition()
    cromer_x.append(position[0])
    cromer_y.append(position[1])

    position = k.averagePosition()
    verlet_x.append(position[0])
    verlet_y.append(position[1])

    position = l.averagePosition()
    rk4_x.append(position[0])
    rk4_y.append(position[1])

for i in data:
    position = i.averagePosition()
    synchro_x.append(position[0])
    synchro_y.append(position[1])

magneticX, magneticY = np.meshgrid(list(range(-3,4,1)), list(range(-5,2,1)))

plt.figure('Euler')
plt.plot(euler_x,euler_y)
plt.scatter(euler_x[-1],euler_y[-1],zorder=10,label='proton bunch')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.axvspan(-0.098, 0.098, color='grey', alpha=0.5 ,zorder=1, label='Electric Field')
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.legend(loc='upper left', framealpha=1)
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Euler-Cromer')
plt.plot(cromer_x,cromer_y)
plt.scatter(cromer_x[-1],cromer_y[-1],zorder=10,label='proton bunch')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.axvspan(-0.098, 0.098, color='grey', alpha=0.5 ,zorder=1, label='Electric Field')
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.legend(loc='upper left', framealpha=1)
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('Velocity Verlet')
plt.plot(verlet_x,verlet_y)
plt.scatter(verlet_x[-1],verlet_y[-1],zorder=10,label='proton bunch')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.axvspan(-0.098, 0.098, color='grey', alpha=0.5 ,zorder=1, label='Electric Field')
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.legend(loc='upper left', framealpha=1)
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.figure('RK4')
plt.plot(rk4_x,rk4_y)
plt.scatter(rk4_x[-1],rk4_y[-1],zorder=10,label='proton bunch')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.axvspan(-0.098, 0.098, color='grey', alpha=0.5 ,zorder=1, label='Electric Field')
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.legend(loc='upper left', framealpha=1)
plt.tick_params(which='both',direction='in',right=True,top=True)

magneticX, magneticY = np.meshgrid(np.linspace(-0.6,0.6,6),np.linspace(-1.2,0.,6))

plt.figure('Synchrocylotron')
plt.plot(synchro_x,synchro_y)
plt.scatter(synchro_x[-1],synchro_y[-1],zorder=10,label='proton bunch')
plt.scatter(magneticX,magneticY,marker=r'$\odot$',s=95,color='black',label='Magnetic Field',zorder=3)
plt.axvspan(-0.05, 0.05, color='grey', alpha=0.5 ,zorder=1, label='Electric Field')
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.legend(loc='upper right', framealpha=1)
plt.tick_params(which='both',direction='in',right=True,top=True)

plt.show()
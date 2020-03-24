import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib

try:
    simulation_data = np.load('cyclotron_data.npz',allow_pickle=True)
except FileNotFoundError:
    import RecordCyclotron
    simulation_data = np.load('cyclotron_data.npz',allow_pickle=True)

time = simulation_data['time']
eulerData = simulation_data['euler']
cromerData = simulation_data['cromer']
verletData = simulation_data['verlet']
RK4Data = simulation_data['rk']

euler_x, euler_y = [], []
cromer_x, cromer_y = [], []
verlet_x, verlet_y = [], []
rk4_x, rk4_y = [], []
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

plt.figure('Euler')
plt.plot(euler_x,euler_y)

plt.figure('Euler-Cromer')
plt.plot(cromer_x,cromer_y)

plt.figure('Velocity Verlet')
plt.plot(verlet_x,verlet_y)

plt.figure('RK4')
plt.plot(rk4_x,rk4_y)
plt.show()
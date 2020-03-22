import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

def compare():
    field = EMField([0,0,0], [0,0,1.6*10**(-5)], [0,0]) 
    protons_1 = ProtonBunch(0.0047,1) ; protons_2 = copy.deepcopy(protons_1)
    protons_3 = copy.deepcopy(protons_1) ; protons_4 = copy.deepcopy(protons_1)
    
    time, deltaT, duration = 0, 10**(-5), 0.0042*10

    inital_bunch_1 = copy.deepcopy(protons_1) ; inital_bunch_2 = copy.deepcopy(protons_2)
    inital_bunch_3 = copy.deepcopy(protons_3) ; inital_bunch_4 = copy.deepcopy(protons_4)

    timeSeries = [0.]
    eulerData = [inital_bunch_1] ; cromerData = [inital_bunch_2]
    verletData = [inital_bunch_3] ; RK4Data = [inital_bunch_4]

    while time <= duration:
        time += deltaT
        timeSeries.append(time)
        field.getAcceleration(protons_1.bunch, time, deltaT)
        field.getAcceleration(protons_2.bunch, time, deltaT)
        field.getAcceleration(protons_3.bunch, time, deltaT)
        field.getAcceleration(protons_4.bunch, time, deltaT)
        protons_1.update(deltaT,field,time,0)
        protons_2.update(deltaT,field,time,1)
        protons_3.update(deltaT,field,time,2)
        protons_4.update(deltaT,field,time,3)
        temp_bunch_1 = copy.deepcopy(protons_1)
        temp_bunch_2 = copy.deepcopy(protons_2)
        temp_bunch_3 = copy.deepcopy(protons_3)
        temp_bunch_4 = copy.deepcopy(protons_4)
        eulerData.append(temp_bunch_1)
        cromerData.append(temp_bunch_2)
        verletData.append(temp_bunch_3)
        RK4Data.append(temp_bunch_4)
    np.savez('methodsData',time=timeSeries,euler=eulerData,cromer=cromerData,verlet=verletData,rk=RK4Data)
    return np.load('methodsData.npz',allow_pickle=True)

try:
    simulation_data = np.load('methodsData.npz', allow_pickle=True)
except FileNotFoundError:
    simulation_data = compare()

time = simulation_data['time']
eulerData = simulation_data['euler']
cromerData = simulation_data['cromer']
verletData = simulation_data['verlet']
RK4Data = simulation_data['rk']

plt.figure('positions')
plt.plot([i.averagePosition()[0] for i in eulerData], [i.averagePosition()[1] for i in eulerData], label='euler')
plt.plot([i.averagePosition()[0] for i in cromerData], [i.averagePosition()[1] for i in cromerData], label='euler cromer')
plt.plot([i.averagePosition()[0] for i in verletData], [i.averagePosition()[1] for i in verletData], label='verlet')
plt.plot([i.averagePosition()[0] for i in RK4Data], [i.averagePosition()[1] for i in RK4Data], label='4th Order RK')
plt.legend()

plt.figure('fractional energies')
plt.plot(time, [i.KineticEnergy()/eulerData[0].KineticEnergy() for i in eulerData], label='euler')
plt.plot(time, [i.KineticEnergy()/cromerData[0].KineticEnergy() for i in cromerData], label='euler cromer')
plt.plot(time, [i.KineticEnergy()/verletData[0].KineticEnergy() for i in verletData], label='verlet')
plt.plot(time, [i.KineticEnergy()/RK4Data[0].KineticEnergy() for i in RK4Data], label='4th Order RK')
plt.legend()
plt.show()
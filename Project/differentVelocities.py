import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([0,0,0], [0,0,0.07], [0,0]) 
protons_1 = ProtonBunch(10**(6),1) ; protons_2 = ProtonBunch(3*10**(6),1)
protons_3 = ProtonBunch(5*10**(6),1) 

time, deltaT, duration = 0, 10**(-9), 10**(-6)*1

timeSeries = [0.]
bunch_1 = [protons_1.averagePosition()] ; bunch_2 = [protons_2.averagePosition()]
bunch_3 = [protons_3.averagePosition()]

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    timeSeries.append(time)
    field.getAcceleration(protons_1.bunch, time, deltaT)
    field.getAcceleration(protons_2.bunch, time, deltaT)
    field.getAcceleration(protons_3.bunch, time, deltaT)

    protons_1.update(deltaT,field,time,2)
    protons_2.update(deltaT,field,time,2)
    protons_3.update(deltaT,field,time,2)

    temp_bunch_1 = copy.deepcopy(protons_1)
    temp_bunch_2 = copy.deepcopy(protons_2)
    temp_bunch_3 = copy.deepcopy(protons_3)

    bunch_1.append(temp_bunch_1.averagePosition())
    bunch_2.append(temp_bunch_2.averagePosition())
    bunch_3.append(temp_bunch_3.averagePosition())

data1_x, data1_y = [], []
data2_x, data2_y = [], []
data3_x, data3_y = [], []
for i,j,k in zip(bunch_1,bunch_2,bunch_3):
    data1_x.append(i[0])
    data1_y.append(i[1])

    data2_x.append(j[0])
    data2_y.append(j[1])

    data3_x.append(k[0])
    data3_y.append(k[1])

plt.figure('DiffVelocities')
plt.plot(data1_x,data1_y,color='blue',label=r'$E_{k,inital} = 1$ MeV')
plt.plot(data2_x,data2_y,color='red',label=r'$E_{k,inital} = 3$ MeV')
plt.plot(data3_x,data3_y,color='green',label=r'$E_{k,inital} = 5$ MeV')
plt.legend(loc='upper right',framealpha=1)
plt.xlabel('x-position [m]')
plt.ylabel('y-position [m]')
plt.tick_params(which='both',direction='in',right=True,top=True)
plt.show()
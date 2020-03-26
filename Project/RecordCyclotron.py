import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib
import copy
from EMField import EMField
from ProtonBunch import ProtonBunch

field = EMField([.1,0,0], [0,0,1.6*10**(-5)], [-0.104,0.104]) 
protons_1 = ProtonBunch(0.0047,5) ; protons_2 = copy.deepcopy(protons_1)
protons_3 = copy.deepcopy(protons_1) ; protons_4 = copy.deepcopy(protons_1)

log.logger.info('Initial average kinetic energy: %s eV' % protons_1.KineticEnergy())
log.logger.info('Initial average momentum: %s kg m/s' % protons_1.momentum())
log.logger.info('Initial average position %s m' % protons_1.averagePosition())
log.logger.info('Initial bunch position spread: %s m' % protons_1.positionSpread())
log.logger.info('Initial bunch energy spread: %s eV' % protons_1.energySpread())
    
time, deltaT, duration = 0, 10**(-5), 0.0041*5

inital_bunch_1 = copy.deepcopy(protons_1) ; inital_bunch_2 = copy.deepcopy(protons_2)
inital_bunch_3 = copy.deepcopy(protons_3) ; inital_bunch_4 = copy.deepcopy(protons_4)

timeSeries = [0.]
eulerData = [inital_bunch_1] ; cromerData = [inital_bunch_2]
verletData = [inital_bunch_3] ; RK4Data = [inital_bunch_4]

log.logger.info('starting simulation')
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
log.logger.info('simulation finished')

log.logger.info('writing lists to file')
np.savez('cyclotron_data',time=timeSeries,euler=eulerData,cromer=cromerData,verlet=verletData,rk=RK4Data)
log.logger.info('file writing complete')
import log
import numpy as np
from Particle import Particle
from ChargedParticle import ChargedParticle
from EMField import EMField
import scipy.constants as const
import matplotlib.pyplot as plt
import copy

EulerProton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [0,3000,0])
EulerCromerProton = ChargedParticle('proton-2', const.m_p, const.e, [0,0,0], [0,3000,0])
field = EMField([0,0,0], [0,0,1.6*10**(-5)])

time, deltaT, duration = 0, 10**(-5), 0.0041

EulerData = []
EulerCromerData = []

log.logger.info('starting simulation')
while time <= duration:
    time += deltaT
    field.getAcceleration(EulerProton)
    field.getAcceleration(EulerCromerProton)
    EulerProton.euler(deltaT)
    EulerCromerProton.eulerCromer(deltaT)
    temp_particle1, temp_particle2 = copy.deepcopy(EulerProton), copy.deepcopy(EulerCromerProton)
    EulerData.append(temp_particle1.position)
    EulerCromerData.append(temp_particle2.position)
log.logger.info('simulation finished')

log.logger.info('creating plots')
Euler_x, Euler_y = [i[0] for i in EulerData], [i[1] for i in EulerData]
Cromer_x, Cromer_y = [i[0] for i in EulerCromerData], [i[1] for i in EulerCromerData]

plt.figure('Constant Magnetic Field')
plt.plot(Euler_x,Euler_y,color='blue',label="Euler Path")
plt.plot(Cromer_x,Cromer_y,color='red',label="Euler Cromer Path")
plt.axvspan(1.5,2.5,alpha=0.5, color='grey')
plt.text(1.7,2.5,'Electric \n Field',fontsize=10,fontweight='bold')
plt.text(0.25,0,'NB: No electric field has been simulated \n       here, the grey section is just POC.',fontsize=10,)
plt.xlim(-1,5)
plt.ylim(-3,3)
plt.xlabel(r'x  position  [$m$]')
plt.ylabel(r'y  position [$m$]')
plt.title('Plot to show the path of a proton in a constant magnetic field')
plt.legend()
plt.show()
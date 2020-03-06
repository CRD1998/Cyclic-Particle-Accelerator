import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
import math
from EMField import EMField
from ChargedParticle import ChargedParticle
import pandas as pd

"""
foo
"""
field = EMField([0,0,0], [0,0,1.6*10**(-5)])
proton = ChargedParticle('proton', const.m_p, const.e, [0,0,0], [3000,0,0])
theoretical_period = 2*const.pi*proton.mass / (proton.charge*field.magneticMag())

def simulation():
    bunch = [proton]

    time, deltaT, duration = 0, 10**(-6), 0.04

    timeSeries = []
    data_csv = []
    log.logger.info('starting simulation')
    while time <= duration:
        time += deltaT
        timeSeries.append(time)
        field.getAcceleration(bunch, time, deltaT)

        temp_bunch = copy.deepcopy(bunch)

        data_csv.append([time,temp_bunch[0].position,temp_bunch[0].velocity,temp_bunch[0].KineticEnergy(),temp_bunch[0].momentum()])
    
    log.logger.info('simulation finished')

    log.logger.info('converting to pandas dataframe')
    df = pd.DataFrame(data_csv, columns=['Time','Position','Velocity','Kinetic Energy','Momentum'])
    df.to_pickle('velocity_data.csv')
    log.logger.info('data pickled')
    
    return pd.read_pickle('velocity_data.csv')

try:
    df = pd.read_pickle('velocity_data.csv')
except FileNotFoundError as err:
    df = simulation()

print(df)
timeSeries = df['Time'].tolist()
velocities = df['Velocity'].tolist()
kineticEnergies = df['Kinetic Energy'].tolist()
inital_speed = np.linalg.norm(velocities[0])
fractional_speeds = [np.linalg.norm(i)/inital_speed for i in velocities]
inital_kinetic = kineticEnergies[0]
fractional_kinetic = [i/inital_kinetic for i in kineticEnergies]
revolutions = np.linspace(0,round(timeSeries[-1]/0.0041),len(fractional_speeds))

plt.figure('Fractional Speeds')
plt.plot(timeSeries,fractional_speeds)
plt.xlabel('Time [s]')
plt.ylabel('Fractional Linear Speed [m/s]')
ax2 = plt.twiny()
ax2.set_ylabel('Revolutions')
ax2.plot(revolutions,fractional_speeds)

plt.figure('Fractional Kinetic Energies')
plt.plot(timeSeries,fractional_kinetic)
plt.xlabel('Time [s]')
plt.ylabel('Fractional Kinetic Energy [J]')

# TODO add fractional momentum as well
#plt.tight_layout() 
plt.show()


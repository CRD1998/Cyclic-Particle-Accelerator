import log
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
import copy
from EMField import EMField
from ChargedParticle import ChargedParticle
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset

"""
This file will plot the fractional linear speed, kinetic energy and momentum against time and the number
of revolutions. A pickled csv file called "velocity_data.csv" is required to plot these graphs. If you
do not have this file the simulation function will run and generate the csv for you.

The theoretical time period is being used here, to see the justification of using the theoreitcal time
period see "analyse_period.csv". This file demonstrates how accurate the simulated time period is.
"""
field = EMField([0,0,0], [0,0,1.6*10**(-5)], [0,0,0]) # join the dees so the magneitc field in constant
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

timeSeries = df['Time'].tolist()
revolutions = np.linspace(0,round(timeSeries[-1]/0.0041),len(timeSeries))

velocities = df['Velocity'].tolist()
kineticEnergies = df['Kinetic Energy'].tolist()
momenta = df['Momentum'].tolist()
inital_speed = np.linalg.norm(velocities[0])
fractional_speeds = [np.linalg.norm(i)/inital_speed for i in velocities]
inital_kinetic = kineticEnergies[0]
fractional_kinetic = [i/inital_kinetic for i in kineticEnergies]
inital_momentum = np.linalg.norm(momenta[0])
fractional_momenta = [np.linalg.norm(i)/inital_momentum for i in momenta]

plt.figure('Fractional Linear Speed')
plt.plot(timeSeries,fractional_speeds)
plt.xlabel('Time [s]')
plt.ylabel('Fractional Linear Speed [m/s]')
ax2 = plt.twiny()
ax2.set_xlabel('Revolutions')
ax2.plot(revolutions,fractional_speeds)

plt.figure('Fractional Kinetic Energies')
plt.plot(timeSeries,fractional_kinetic)
plt.xlabel('Time [s]')
plt.ylabel('Fractional Kinetic Energy [J]')
ax2 = plt.twiny()
ax2.set_xlabel('Revolutions')
ax2.plot(revolutions,fractional_kinetic)

plt.figure('Fractional Momentum')
plt.plot(timeSeries,fractional_momenta)
plt.xlabel('Time [s]')
plt.ylabel('Fractional Momentum [J]')
ax2 = plt.twiny()
ax2.set_xlabel('Revolutions')
ax2.plot(revolutions,fractional_momenta)

plt.show()
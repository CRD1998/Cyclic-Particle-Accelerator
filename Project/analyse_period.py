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
This file is an independent file that tests the simulated time period of a proton in a constant magnetic 
field against the expected time period given by:
T = 2*pi*m / (qB)
If you do not have the csv containing approximately 100 revolutions, the simulation will run and generate
the csv file for you. It will then print the mean time period and its error (one standard deviation) and the
time period predicted by the equation above. Finally it will print, the ratio of the two time periods; the 
simulated time period divided by the theoretical one.
"""
field = EMField([0,0,0], [0,0,1.6*10**(-5)])
proton = ChargedParticle('proton', const.m_p, const.e, [0,0,0], [3000,0,0])
theoretical_period = 2*const.pi*proton.mass / (proton.charge*field.magneticMag())

def simulation():
    bunch = [proton]

    time, deltaT, duration = 0, 10**(-5), 0.4

    timeSeries = []
    data_csv = []
    revolution = 1
    data_csv.append([revolution,time,0,[0,0,0]]) # ensure data always has one element to avoid length error later
    log.logger.info('starting simulation')
    while time <= duration:
        time += deltaT
        timeSeries.append(time)
        previous_proton = copy.deepcopy(bunch)
        field.getAcceleration(bunch, time, deltaT)

        temp_bunch = copy.deepcopy(bunch)

        if previous_proton[0].position[0] < 0 and temp_bunch[0].position[0] >= 0:
            period = time - (data_csv[-1])[1]
            revolution += 1
            log.logger.info('computing revolution')
            data_csv.append([revolution,time,period,temp_bunch[0].position])
    
    log.logger.info('simulation finished')

    log.logger.info('converting to pandas dataframe')
    df = pd.DataFrame(data_csv, columns=['Revolution','Time','Measured Period','Position'])
    df.to_pickle('period_data.csv')
    log.logger.info('data pickled')
    
    return pd.read_pickle('period_data.csv')

try:
    df = pd.read_pickle('period_data.csv')
except FileNotFoundError as err:
    df = simulation()

print(df)
df.drop([df.index[0]], inplace=True) 
average_period = df['Measured Period'].mean() # seconds
error = df['Measured Period'].std() # seconds
print('time period from simulation: {0} ± {1} s'.format(round(average_period,6),round(error,6)))
print('time period from theory: {0} s'.format(round(theoretical_period,6)))
print('raito of simulated time period to theoretical: {0}'.format(round((average_period/theoretical_period),6)))
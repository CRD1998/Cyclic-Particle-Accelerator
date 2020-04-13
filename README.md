# phys389-2020-project-CRD1998
phys389-2020-project-CRD1998 created by GitHub Classroom

The simulation attached models a proton bunch in a synchrocyclotron. There are several files required to run this simulation, which have 
all been detailed below.

Simulation Files
----------------

1. Particle.py - Conatins a class also called Particle, generates an object to represent a subatomic paritcle, such as a neutron. This 
                 class also contains methods returning the Lorentz factor, momentum, kinetic energy and methods to update an instances
                 position using one of the following methods: Euler, Euler-Cromer, velocity Verlet, fourth order Runge-Kutta.
2. ChargedParticle.py - Inherits from the class Particle, generates an object to represent a charged subatomic particle, such as a proton.
3. EMField.py - Contains a class called EMField, this generates an object to represent to an electromagnetic field in a (synchro)cyclotron.
                It consists of a constant magnetic field and a time-varying electric field in the form Acos(wt+phi) bounded by two x
                ordinates.
4. Bunch.py - An abstract base class that provides the blueprint for an object representing a bunch of ChargedParticle objects, for example
              a bunch of protons. It samples initial positions and velocities for these particles from a Gaussian and can return the Lorentz
              factor of the bunch, average position and average velocity. It also contains a method to reduce the time step used in the
              update methods if any of the particles in the bunch are near or in the electric field, to ensure the electric field is not
              stepped over.
5. ProtonBunch.py - This files contains a class called ProtonBunch, which inherits from the Bunch ABC, it represents a bunch of protons.
6. RecordCylotron/Synchrocyclotron.py - This is the main simulation file and should be the only one a user interacts with. This is where 
                                        the particle bunches and EM fields are instantitated, as well as where the time step, duration of
                                        the simulation and theupdate to be used are declared. Once the simulation has finished, a list of 
                                        the time values is written to an npz file alongside a list containing a copy of the ProtonBunch 
                                        object at every time value.

Analysis Files
--------------

1. analyse_methods.py - First, this file will look for a npz file called "methods_data.npz" which contains data from simulating a proton
                        in a constant magnetic field for approximately 100 revolutions. It will then plot the fractional kinetic energy 
                        and momentum of this proton against time. If the "methods_data.npz" is not found, the data file from RecordCyclotron.py
                        will searched for, if this is also not found it will be generated and the same plots will be produced. In the case
                        the data file from RecordCyclotron.py needs to be generated, ensure you remove the electric field entirely by 
                        setting its field strength and width to zero so only a constant magnetic field is present. By only have a magnetic 
                        field present, both energy and momentum should be constant throughout the simulation.
2. analyse_period.py - This script will first look for a pickled csv named "period_data.csv" which contains a pandas dataframe detailing
                       how many revolutions were measured when simulating a proton in a constant magnetic field. It then averages the 
                       measured time period and take ones standard deviation as its error, which it comapres to the theoretical time period.
                       If you do not have this csv, it will be generated for you.
3. analyse_phase.py - This script simulates a bunch of protons in a cyclotron for approximately 10 revolutions nine times, each iteration
                      has a different phase shift applied to the electric field. It then plots the spread in the y-position and the spread
                      in kineitc energy against time for each phase shift. It first looks for a file containing this data called "phase_data.npz"
                      if you do not have this file it will be generated for you.
4. analyse_synchro_phase.py - This script simulates a bunch of protons in a synchrocyclotron for approximately 10 revolutions nine times, 
                              each iteration has a different phase shift applied to the electric field. It then plots the spread in the 
                              y-position and the spread in kineitc energy against time for each phase shift, it also fits a linear curve 
                              to kinetic energy spread data, so the effects of the phase shift are more obvious. It first looks for a file 
                              containing this data called "phase_synchro_data.npz" if you do not have this file it will be generated for you.
5. analyse_timestep.py - This script simulates a bunch of protons in a cyclotron for approximately 100 revolutions three times, each iteration
                         has a different time step being used in the integrator. It then plots the fractional kinetic energy and momentum against
                         time. It first looks for a file containing this data called "timestep_data.npz"if you do not have this file it will be 
                         generated for you.

Test Files
----------

1. test_ChargedParticles.py - This contains several unit tests. Ensuring the methods in this class, including the ones it inherits are
                              working as expected.
2. test_EMField.py - This contains several unit tests. Ensuring the methods in this class are working as expected. 
3. test_ProtonBunch.py - This contains several unit tests. Ensuring the methods in this class, including the ones it inherits are
                         working as expected.

Notes
-----

A file named "log.py" is also included in this project, this writes to a file named "PHYS389.log" and can be used to check on the progress
of longer simulations as well as follow the execution order of processes in the simulation. The data files required to run all of the 
analysis files have been inlcuded. If you wish to generate your own, I recommend reducing the simulation's duration and reducing the number
of particles in the bunch to prevent long run times.

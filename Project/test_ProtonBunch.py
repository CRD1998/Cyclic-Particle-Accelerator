import pytest
import scipy.constants as const
import numpy as np
from unittest.mock import patch
from ProtonBunch import ProtonBunch
from ChargedParticle import ChargedParticle

proton_1 = ChargedParticle('proton-1', const.m_p, const.e, [50,-20,10], [1000,-2000,3000])
proton_2 = ChargedParticle('proton-2', const.m_p, const.e, [-75,-30,80], [4000,5000,-6000])
bunch_of_protons = ProtonBunch(3,2) # initialise a ProtonBunch
bunch_of_protons.bunch = [proton_1,proton_2] # overwrite the bunch attribute to replace pseudo-random ChargedParticle objects

def test_averagePosition():
    calculated_position = np.array([-12.5,-25,45],dtype=float)
    assert np.array_equal(calculated_position, bunch_of_protons.averagePosition())

def test_averageVelocity():
    calculated_velocity = np.array([2500,1500,-1500],dtype=float)
    assert np.array_equal(calculated_velocity, bunch_of_protons.averageVelocity())

def test_momentum():
    """
    This will test that the momentum function in the Bunch ABC is correctly returning the average
    kinetic energy for a particle in the bunch. If method is called and True is parsed into it, then
    then this will also test that the method is correctly returning the total momentum of the bunch.
    """
    calculated_average_momentum = const.m_p * np.array([2500,1500,-1500],dtype=float)
    calculated_total_momentum = const.m_p * np.array([5000,3000,-3000],dtype=float)

    assert np.allclose(calculated_average_momentum,bunch_of_protons.momentum(),rtol=10**(-9))
    assert np.allclose(calculated_total_momentum,bunch_of_protons.momentum(True),rtol=10**(-9))

def test_KineticEnergy():
    """
    This will test that the KineticEnergy function in the Bunch ABC is correctly returning the 
    average kinetic energy for a particle in the bunch. If method is called and True is parsed into 
    it, then then this will also test that the method is correctly returning the total kinetic 
    energy of the bunch.
    """
    calculated_average_kinetic = 0.237502831 # eV
    calculated_total_kinetic = 0.475005663 # eV

    assert calculated_average_kinetic == pytest.approx(bunch_of_protons.KineticEnergy())
    assert calculated_total_kinetic == pytest.approx(bunch_of_protons.KineticEnergy(True))

@patch.object(ProtonBunch, 'distributeEnergies')
def test_assignVelocities(mock_energies):
    """
    This will test that the assignVelocities method correctly calculates a linear speed given a
    kinetic energy. It also tests that the function returns a list of velocity vectors where the
    calculated speeds are the x-components in these velocity vectors.
    """
    mock_energies.return_value=np.array([10**8,0.75]) # return these values instead of values sampled from normal distribution
    velocity_1, velocity_2 = [128369776.9,0,0], [11986.4508,0,0]
    calculated_velocities = [velocity_1, velocity_2]
    assigned_velocities = bunch_of_protons.assignVelocities()
    assert np.allclose(calculated_velocities, assigned_velocities, rtol=10**(-4)) # FIXME Why only accurate to 10^(-4)?

def test_spread():
    """
    This will test that the spread method (numpy standard deviation) is correctly returning the variance
    of all the particles' positions in a bunch.
    """
    calculated_spread = np.array([62.5,5,35],dtype=float)
    assert np.array_equal(calculated_spread, bunch_of_protons.spread())

@patch.object(ProtonBunch, 'assignVelocities')
@patch.object(ProtonBunch, 'assignPositions')
def test_createBunch(mock_positions, mock_velocities):
    mock_positions.return_value = [[1,2,3], [4,5,6]]
    mock_velocities.return_value = [[10,20,30], [40,50,60]]
    particle_1 = ChargedParticle('proton-1', const.m_p, const.e, [1,2,3], [10,20,30])
    particle_2 = ChargedParticle('proton-2', const.m_p, const.e, [4,5,6], [40,50,60])
    expected_bunch = [particle_1, particle_2]
    actual_bunch = bunch_of_protons.createBunch()
    assert np.array_equal(expected_bunch, actual_bunch)
    
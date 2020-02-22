import pytest
from Particle import Particle
import scipy.constants as const
import numpy as np

def test_exceedingC():
    """
    This test will check that if the user tries to create a particle whose speed exceeds the speed
    of light in a vacuum, a ValueError is raised.
    """
    with pytest.raises(ValueError):
        Particle('proton', const.m_p, [0,0,0], [300000000,0,0])

def test_gamma():
    """
    This test checks that the gamma function correctly calculates the Lorenzt factor. Using the 
    approximation function accounts for the floating point error. 
    """
    calculated_value = 1.342384701 # calculated by hand using c = 299792458 [ms^(-1)]
    proton = Particle('proton', const.m_p, [0,0,0], [200000000,0,0])
    assert calculated_value == pytest.approx(proton.gamma())

def test_Euler():
    """
    This test checks that the Euler method is correctly updating a particle's position and
    velocity.
    """
    proton = Particle('proton', const.m_p, [0,0,0], [-100,200,300], [10,20,-30])
    deltaT = 0.5
    calculated_position = np.array([-50,100,150], dtype=float)
    calculated_velocity = np.array([-95,210,285], dtype=float)
    proton.euler(deltaT)
    assert np.array_equal(calculated_position,proton.position)
    assert np.array_equal(calculated_velocity,proton.velocity)

def test_EulerCromer():
    """
    This test checks that the Euler Cromer method is correctly updating a particle's position and
    velocity.
    """
    proton = Particle('proton', const.m_p, [0,0,0], [-100,200,300], [10,20,-30])
    deltaT = 0.5
    calculated_velocity = np.array([-95,210,285], dtype=float)
    calculated_position = np.array([-47.5,105,142.5], dtype=float)
    proton.eulerCromer(deltaT)
    assert np.array_equal(calculated_velocity,proton.velocity)
    assert np.array_equal(calculated_position,proton.position) 
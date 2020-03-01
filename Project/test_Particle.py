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

@pytest.mark.parametrize('test_input,expected',[(Particle('proton', const.m_p, [0,0,0], [3000,0,0]),7.526911023*10**(-21)),
                        (Particle('proton', const.m_p, [0,0,0], [200000000,0,0]),5.146992568*10**(-11))])
def test_KineticEnergy(test_input,expected):
    """
    This test checks that the Kinetic Energy method correctly returns the kinetic energy 
    of a particle in both the relativistic and non-relativistic case.
    """    
    assert expected == pytest.approx(test_input.KineticEnergy())

@pytest.mark.parametrize('test_input,expected',[(Particle('proton', const.m_p, [0,0,0], [3000,0,0]),const.m_p*np.array([3000,0,0],dtype=float)),
                        (Particle('proton', const.m_p, [0,0,0], [200000000,0,0]),const.m_p*1.342384701*np.array([200000000,0,0],dtype=float))])
def test_momentum(test_input,expected):
    """
    This test checks that the momentum method correctly returns the momentum 
    of a particle in both the relativistic and non-relativistic case.
    """    
    assert np.allclose(test_input.momentum(),expected,rtol=10**(-10))

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
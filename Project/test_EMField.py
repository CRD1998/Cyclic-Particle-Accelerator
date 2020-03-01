import pytest
from ChargedParticle import ChargedParticle
from EMField import EMField
import scipy.constants as const
import numpy as np

def test_frequency():
    """
    This function tests that the cyclotron frequency is correctly being calculated.
    """
    proton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [2000,-4000,6000])
    field = EMField([-0.3,0.1,0.2], [6*10**(-5),-7*10**(-5),8*10**(-5)])
    calculated_frequency = 11692.45605
    assert calculated_frequency == pytest.approx(field.frequency(proton))


@pytest.mark.parametrize('test_input,expected',
                        [((ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [2000,0,0]),EMField([0.1,0,0])), [0.1,0,0]),
                        ((ChargedParticle('proton-1', const.m_p, const.e, [10,0,0], [2000,0,0]),EMField([0.1,0,0])), [0,0,0]),
                        ((ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [2000,-4000,6000]),EMField([-0.3,0.1,0.2], [6*10**(-5),-7*10**(-5),8*10**(-5)])), [-0.2,0.3,0.3])])
def test_getAcceleration(test_input,expected):
    """
    This function will test the getAcceleration method in the EMField class. It will check that in
    the presence of both an electric and magnetic field, the acceleration is correctly calculated. 
    It will also check that the electric field is being correctly constrained. The three test cases
    are:
        1) An electric field exists within the specified coordinates
        2) The electric field is zero outside of these coordinates.
        3) Given both an electric and magnetic field, the acceleration due to the Lorentz force
           is correctly being calculated.
    """
    proton, field = test_input
    field.getAcceleration(proton,0) # set time to 0 so the electric field is equal to its magnitude 
    calculated_lorentz = np.array(expected,dtype=float) # get Lorentz force
    calculated_lorentz *= const.physical_constants['proton charge to mass quotient'][0] # convert to acceleration
    assert np.allclose(calculated_lorentz, proton.acceleration, rtol=10**(-9))
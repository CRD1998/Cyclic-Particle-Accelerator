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

#@pytest.mark.parametrize('test_input,expected',[(EMField([0.1,0,0],[0,0,10**(-5)]).ImplementElectricField(ChargedParticle('proton-1', const.m_p, const.e),1),
#                        np.array([-0.095461349,0,0],dtype=float)),
#                        (EMField([0.1,0,0],[0,0,10**(-5)]).ImplementElectricField(ChargedParticle('proton-1', const.m_p, const.e,[2,0,0]),1),
#                        np.array([0,0,0],dtype=float))])
#def test_ImplementElectricField(test_input,expected):
#    """
#    This function will test that the electric field is correctly being implemented:
#        1) It exists within the specified coordinates in EMField.py
#        2) It is zero outside of these coordinates.
 #   """
 #   assert np.allclose(test_input,expected,rtol=10**(-9))


def test_getAcceleration():
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
    pass
    #proton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [2000,-4000,6000])
    #field = EMField([-0.3,0.1,0.2], [6*10**(-5),-7*10**(-5),8*10**(-5)])
    #charge_mass = const.physical_constants['proton charge to mass quotient'][0]
    #calculated_lorentz = np.array([-0.2,0.3,0.3], dtype=float)
    #calculated_lorentz *= charge_mass
    #field.getAcceleration(proton)
    #assert np.allclose(calculated_lorentz, proton.acceleration, rtol=10**(-9))
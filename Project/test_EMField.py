import pytest
from ChargedParticle import ChargedParticle
from EMField import EMField
import scipy.constants as const
import numpy as np

def test_getAcceleration():
    proton = ChargedParticle('proton-1', const.m_p, const.e, [0,0,0], [2000,-4000,6000])
    field = EMField([-0.3,0.1,0.2], [6*10**(-5),-7*10**(-5),8*10**(-5)])
    charge_mass = const.physical_constants['proton charge to mass quotient'][0]
    calculated_lorentz = np.array([-0.2,0.3,0.3], dtype=float)
    calculated_lorentz *= charge_mass
    field.getAcceleration(proton)
    assert np.allclose(calculated_lorentz, proton.acceleration, rtol=10**(-9))
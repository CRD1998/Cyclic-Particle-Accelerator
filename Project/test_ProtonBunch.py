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
    mock_energies.return_value=np.array([10**8,0.75])
    velocity_1, velocity_2 = [128369776.9,0,0], [11986.4508,0,0]
    calculated_velocities = [velocity_1, velocity_2]
    assigned_velocities = bunch_of_protons.assignVelocities()
    assert np.allclose(calculated_velocities, assigned_velocities, rtol=10**(-4))
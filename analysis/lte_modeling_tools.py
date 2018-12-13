import numpy as np
from astropy import constants
from astropy.modeling.models import custom_model
from astropy import units as u

def nupper_of_kkms(kkms, freq, Aul, degeneracies, replace_bad=None):
    """ Derived directly from pyspeckit eqns..."""

    if replace_bad:
        neg = kkms <= 0
        kkms[neg] = replace_bad

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)
    #nline = 1.95e3 * freq**2 / Aul * kkms
    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2
    # term2 = np.exp(-constants.h*freq/(constants.k_B*Tex)) -1
    # term2 -> kt / hnu
    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)
    return (nline * Khz / degeneracies).to(u.cm**-2)

def kkms_of_nupper(nupper, freq, Aul, degeneracies):
    """
    Convert the column density in the upper state of a line ``nupper`` to the
    integrated intensity in brightness units (K km / s).
    """

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    nupper = u.Quantity(nupper, u.cm**-2)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    Khz = (nupper / nline * degeneracies)

    kkms = (Khz / (freq/constants.c)).to(u.K * u.km/u.s)

    return kkms


def rovib_lte_model_generator(vibenergies, rotenergies):

    @custom_model
    def model(jstate, vstate, logcolumn=np.log(1e13), rottem=100, vibtem=2000):
        elower_vib = np.array([vibenergies[int(v)] for v in vstate])
        eupper_j = np.array([rotenergies[ju] for ju in jstate])
        #result = -1/rottem * (eupper - elower_vib) + column - 1/vibtem * eupper

        # these are the populations of states in the v=0 state at a given
        # J-level.
        result_v0 = -1/rottem * eupper_j + logcolumn

        # Then, for each of these, we determine the levels in the vibrationally
        # excited state by adding e^-hnu/kt, where t is a different temperature
        # (t_v) and nu is now just the nu provided by the vibrations
        result = result_v0 - 1/vibtem * elower_vib

        return result

    return model()

def simple_lte_model_generator():

    @custom_model
    def model(eupper, logcolumn=np.log(1e13), tem=100):
        """
        Calculate the quantity N_u/g_u as a function of E_u in Kelvin

        The 'logcolumn' quantity is N_tot / Q_tot

        Temperature is the excitation temperature
        """

        result = -1/tem * eupper + logcolumn

        return result

    return model()



if __name__ == '__main__':
    # round-trip test
    kkms = 100*u.K*u.km/u.s
    freq = 100*u.GHz
    Aul = 1*u.s**-1
    degeneracies = 1
    nupper = nupper_of_kkms(kkms, freq, Aul, degeneracies)
    kkms2 = kkms_of_nupper(nupper, freq, Aul, degeneracies)
    np.testing.assert_almost_equal(kkms2.value, kkms.value)

import numpy as np
from astropy import constants
from astropy.modeling.models import custom_model
from astropy import units as u

def nupper_of_kkms(kkms, freq, Aul, replace_bad=None):
    """
    Mangum & Shirley 2015 eqn 82 gives, for the optically thin, Rayleigh-Jeans,
    negligible background approximation:

        Ntot = (3 k) / (8 pi^3 nu S mu^2 R_i)   (Q/g) exp(E_u/k Tex) integ(T_R/f dv)

    Eqn 31:

        Ntot/Nu = Q_rot / gu exp(E_u/k Tex)

        -> Ntot = Nu Q_rot / gu exp(E_u/k Tex)
        -> Nu = N_tot g / Qrot exp(-E_u / k Tex)

    To get Nu of an observed line, then:

        Nu (Q_rot / gu) exp(E_u/k Tex) = (3 k) / (8 pi^3 nu S mu^2 R_i)   (Q/g) exp(E_u/k Tex) integ(T_R/f dv)

    This term cancels:
        (Q_rot / gu) exp(E_u/k Tex)

    Leaving:

        Nu = (3 k) / (8 pi^3 nu S mu^2 R_i)   integ(T_R/f dv)

    integ(T_R/f dv) is the optically thin integrated intensity in K km/s
    dnu/nu = dv/c [doppler eqn], so to get integ(T_R dnu), sub in dv = c/nu dnu

        Nu = (3 k c) / (8 pi^3 nu^2  S mu^2 R_i)   integ(T_R/f dnu)


    We then need to deal with the S mu^2 R_i term.  We assume R_i = 1, since we
    are not measuring any hyperfine transitions (R_i is the hyperfine
    degeneracy; eqn 75)
    Equation 11:

        A_ul = 64 pi^4 nu^3 / (3 h c^3) |mu_ul|^2

    Equation 62:

        |mu_ul|^2 = S mu^2

        -> S mu^2 = (3 h c^3 Aul) / (64 pi^4 nu^3)

    Plugging that in gives

        Nu = (3 k c) / (8 pi^3 nu^2  ((3 h c^3 Aul) / (64 pi^4 nu^3)))   integ(T_R/f dnu)
           = (3 k c 64 pi^4 nu^3) / (8 pi^3 nu^2 3 h c^3 Aul)            integ(T_R/f dnu)
           = (8 pi nu k / (Aul c^2 h)) integ(T_R/f dnu)

    which is the equation implemented below.  We could also have left this in
    dv units by substituting du = nu/c dv:

           = (8 pi nu^2 k / (Aul c^3 h)) integ(T_R/f dv)

    """

    if replace_bad:
        neg = kkms <= 0
        kkms[neg] = replace_bad

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)

    return (nline * Khz).to(u.cm**-2)

def kkms_of_nupper(nupper, freq, Aul):
    """
    Convert the column density in the upper state of a line ``nupper`` to the
    integrated intensity in brightness units (K km / s).

    Inversion of nupper_of_kkms above.
    """

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    nupper = u.Quantity(nupper, u.cm**-2)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    Khz = (nupper / nline)

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

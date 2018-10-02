import numpy as np
from astropy import units as u
from astropy import constants


def collision_rate_NaCl(Jupper, temperature, Aul):

    # rotational constant (NIST):
    Be = 0.2180630*u.cm**-1
    # Hebert '68 via Barton:
    dipole_moment = 8.9721*u.D
    # equilibrium bond length r_e (Barton):
    #re = 2.360796042 * u.AA


    temperature_K = temperature.to(u.K).value
    # eqns from Dickinson 1975
    #B = 1.24e8 / (re/u.AA * Be/(u.cm**-1) * Jupper)**2

    # C must have units of inverse energy
    C = (1.93e4 / (dipole_moment/u.D * Be/(u.cm**-1) * Jupper) *
         np.exp(-1.18/(dipole_moment/u.D))) * u.eV**-1

    # beta is defined between 2.18 and 2.19
    #beta = 11600 / temperature
    beta = (constants.k_B**-1 / u.K).to(u.eV**-1)

    # deltaE is the 'threshold energy' from 2.19
    deltaE = (2.48e-4 * Be/u.cm**-1 * (Jupper+1)) * u.eV

    collrate = (1.44e-6 / temperature_K**0.5 * Aul * np.exp(-beta * deltaE) *
                np.log(C * deltaE + C / beta *
                       np.exp(-0.577/(1+2*beta*deltaE))) * u.cm**3)

    return collrate


if __name__ == "__main__":

    from salt_tables import NaCl
    from astropy.table import Column

    collrates100 = collision_rate_NaCl(NaCl['Ju'],
                                       100*u.K,
                                       NaCl['Aij'].quantity,
                                      )
    NaCl.remove_column('alpha_100K')
    NaCl.add_column(Column(name='alpha_100K', data=collrates100))

    collrates500 = collision_rate_NaCl(NaCl['Ju'],
                                       500*u.K,
                                       NaCl['Aij'].quantity,
                                      )
    NaCl.remove_column('alpha_500K')
    NaCl.add_column(Column(name='alpha_500K', data=collrates500))

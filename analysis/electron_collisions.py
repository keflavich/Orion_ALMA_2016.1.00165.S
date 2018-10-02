import numpy as np
from astropy import units as u
from astropy import constants


def collision_rate_NaCl(Jupper, temperature):

    # rotational constant (NIST):
    # (I assume this is B0 in Dickinson 1975)
    Be = 0.2180630*u.cm**-1

    # Hebert '68 via Barton:
    dipole_moment = 8.9721*u.D

    # equilibrium bond length r_e (Barton):
    # (not used, so commented out)
    #re = 2.360796042 * u.AA


    # eqn from Dickinson 1975 table 1
    # (B is not used, so it is commented out)
    #B = 1.24e8 / (re/u.AA * Be/(u.cm**-1) * Jupper)**2

    # eqn from Dickinson 1975 table 1
    # C must have units of inverse energy by inference, however I don't know
    # for sure that its representation in Table 1 has units of eV^-1, I just
    # assume it
    C = (1.93e4 / (dipole_moment/u.D * Be/(u.cm**-1) * Jupper) *
         np.exp(-1.18/(dipole_moment/u.D)**3)) * u.eV**-1

    # beta is defined between 2.18 and 2.19 in Dickinson 1975:
    # it is the inverse Boltzmann constant in eV
    #beta = 11600 / temperature
    beta = (constants.k_B**-1 / u.K).to(u.eV**-1)

    # deltaE is the 'threshold energy' from 2.19 in Dickinson 1975
    deltaE = (2.48e-4 * Be/u.cm**-1 * (Jupper+1)) * u.eV

    # Equation 2.17 in Dickinson 1975
    A = 2.470 * (dipole_moment/u.D)**2 * (Jupper / (2*Jupper + 1))

    # finally, equation 2.23 in Dickinson 1975
    collrate = (1.44e-6 / (temperature/u.K)**0.5 * A * np.exp(-beta * deltaE) *
                np.log(C * deltaE + C / beta *
                       np.exp(-0.577/(1+2*beta*deltaE))) * u.cm**3)

    return collrate


if __name__ == "__main__":

    from salt_tables import NaCl
    from astropy.table import Column

    collrates100 = collision_rate_NaCl(NaCl['Ju'],
                                       100*u.K,
                                      )
    NaCl.remove_column('alpha_100K')
    NaCl.add_column(Column(name='alpha_100K', data=collrates100))

    collrates500 = collision_rate_NaCl(NaCl['Ju'],
                                       500*u.K,
                                      )
    NaCl.remove_column('alpha_500K')
    NaCl.add_column(Column(name='alpha_500K', data=collrates500))

import numpy as np
from astropy import units as u
from astropy import constants

# rotational constant (NIST):
# (I assume this is B0 in Dickinson 1975)
NaCl_Be = 0.2180630*u.cm**-1

# Hebert '68 via Barton:
NaCl_dipole_moment = 8.9721*u.D

# NIST
KCl_Be = 0.1286347 * u.cm**-1
KCl_dipole_moment = 10.239 * u.D

def collision_rate(Jupper, temperature, Be=NaCl_Be,
                   dipole_moment=NaCl_dipole_moment):

    Jlower = Jupper - 1


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
    beta = (constants.k_B**-1 / temperature).to(u.eV**-1)

    # deltaE is the 'threshold energy' from 2.19 in Dickinson 1975
    deltaE = (2.48e-4 * Be/u.cm**-1 * (Jlower + 1)) * u.eV

    # Equation 2.17 in Dickinson 1975
    A = 2.470 * (dipole_moment/u.D)**2 * (Jupper / (2*Jlower + 1))

    # finally, equation 2.23 in Dickinson 1975
    collrate = (1.44e-6 / (temperature/u.K)**0.5 * A * np.exp(-beta * deltaE) *
                np.log(C * deltaE + C / beta *
                       np.exp(-0.577/(1+2*beta*deltaE))) * u.cm**3 * u.s**-1)

    return collrate


if __name__ == "__main__":

    from salt_tables import NaCl, KCl
    from astropy.table import Column

    for ii in range(1, 10):
        print(ii, collision_rate(ii, 100*u.K))

    collrates100 = collision_rate(NaCl['Ju'],
                                  100*u.K,
                                  )
    if 'alpha_100K' in NaCl.columns:
        NaCl.remove_column('alpha_100K')
    NaCl.add_column(Column(name='alpha_100K', data=collrates100))

    collrates500 = collision_rate(NaCl['Ju'],
                                  500*u.K,
                                  )
    if 'alpha_500K' in NaCl.columns:
        NaCl.remove_column('alpha_500K')
    NaCl.add_column(Column(name='alpha_500K', data=collrates500))

    for ii in range(1, 10):
        print(ii, collision_rate(ii, 100*u.K, Be=KCl_Be,
                                 dipole_moment=KCl_dipole_moment))

    collrates100 = collision_rate(KCl['Ju'],
                                  100*u.K, Be=KCl_Be,
                                  dipole_moment=KCl_dipole_moment,
                                  )
    if 'alpha_100K' in KCl.columns:
        KCl.remove_column('alpha_100K')
    KCl.add_column(Column(name='alpha_100K', data=collrates100))

    collrates500 = collision_rate(KCl['Ju'],
                                  500*u.K, Be=KCl_Be,
                                  dipole_moment=KCl_dipole_moment,
                                  )
    if 'alpha_500K' in KCl.columns:
        KCl.remove_column('alpha_500K')
    KCl.add_column(Column(name='alpha_500K', data=collrates500))

import numpy as np
from astropy import units as u
from astropy import constants

# ionization potentials
# https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)

e_Na = 5.13908*u.eV
e_K = 4.34066*u.eV
e_Al = 5.98577*u.eV


def ion_density(temperature, ionization_energy, density, degeneracy_upper=3,
                degeneracy_lower=1):
    # Saha eqn

    # de Broglie wavelength
    lam = (constants.h**2 / (2*np.pi*constants.m_e*constants.k_B * temperature))**0.5
    saha = 2 / lam**3 * np.exp(-ionization_energy/(constants.k_B*temperature))

    # n_ion + n_neutral = n_total
    # n_ion * n_electron / (n_total - n_ion) = saha
    # for n_ion=n_electron...
    # n_ion**2 / (n_total - n_ion) = saha
    # n_ion**2 + n_ion*saha - n_total * saha = 0
    # n_ion = (-saha +/- (saha**2 - 4*1*(n_total*saha))**0.5)/2
    n_ion = (-saha + (saha**2 + 4*1*(density*saha))**0.5)/2
    return (n_ion/density).decompose()

if __name__ == '__main__':

    import pylab as pl

    pl.figure(1).clf()

    tems = np.linspace(50, 5000, 1000)*u.K

    for ii, density in ((1,1e4), (2,1e6), (3,1e8)):
        pl.subplot(3,1,ii)
        density=density*u.cm**-3
        pl.plot(tems, ion_density(tems, e_Na, density=density), label='Na')
        pl.plot(tems, ion_density(tems, e_K, density=density), label='K')
        pl.plot(tems, ion_density(tems, e_Al, density=density), label='Al')
        pl.xlabel("Temperature")
        pl.title("Density = {0:0.0e}".format(density))

        print()
        print(density)
        print("Na", tems[np.argmin(np.abs(0.5-ion_density(tems, e_Na, density=density)))])
        print("K", tems[np.argmin(np.abs(0.5-ion_density(tems, e_K, density=density)))])
        print("Al", tems[np.argmin(np.abs(0.5-ion_density(tems, e_Al, density=density)))])

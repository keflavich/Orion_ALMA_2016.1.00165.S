from astropy import constants
from astropy import units as u
import numpy as np

# Some predicted values from the Chiang & Golreich 1997 model:

# eqns from Chiang & Golreich 1997, eqn 4:
luminosity = 1e4 * u.L_sun # nocite
star_temperature = 4000 * u.K # Testi+ 2010
radius = ((constants.sigma_sb * 4 * np.pi * (star_temperature)**4 / (luminosity))**(-1/2)).to(u.R_sun)
disk_inner_rad = 35*u.au
disk_outer_rad = 60*u.au
outertem = ((2/3/np.pi)**(1/4) * (radius/(disk_outer_rad))**(3/4.) * star_temperature).decompose()
innertem = ((2/3/np.pi)**(1/4) * (radius/(disk_inner_rad))**(3/4.) * star_temperature).decompose()

print("Disk inner, outer radius: {0}, {1}".format(disk_inner_rad, disk_outer_rad))
print("Flat Disk inner, outer temperature: {0}, {1}".format(innertem, outertem))


# equilibrium temperatures
# (these are for a radiating surface that faces the central star at a fixed
# distance; technically the above calculation is just a different equilibrium)
innertem_equilibrium = ((constants.sigma_sb * 4 * np.pi * (disk_inner_rad)**2 / (luminosity))**(-1/4)).decompose()
outertem_equilibrium = ((constants.sigma_sb * 4 * np.pi * (disk_outer_rad)**2 / (luminosity))**(-1/4)).decompose()

print("Face-on equilibrium inner, outer temperature: {0}, {1}".format(innertem_equilibrium, outertem_equilibrium))

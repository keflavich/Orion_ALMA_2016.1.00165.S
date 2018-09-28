import numpy as np
import paths
import parse_exomol_files
from astropy import units as u
from astropy import constants
from astropy.table import Table, Column

if 'nacl_exomol' not in locals():
    nacl_exomol = parse_exomol_files.ExoMol('NaCl', '23Na-35Cl')

nacltbl = nacl_exomol.transitions

gas_temperature = 100*u.K

stellar_temperature = 4000*u.K
stellar_luminosity = 2e4*u.L_sun
stellar_radius = ((stellar_luminosity / (constants.sigma_sb *
                                         stellar_temperature**4) /
                   (4*np.pi))**0.5).to(u.R_sun)
voldens = 1*u.cm**-3 # assume n(h) ~ 10^10, X(NaCl) ~ 10^-10

pf = Table.read(paths.salty('23Na-37Cl__Barton.pf'), format='ascii')

# partition function: look up from table
Zstar = pf['col2'][np.argmin(np.abs(pf['col1'] - stellar_temperature.value))]
Zgas = pf['col2'][np.argmin(np.abs(pf['col1'] - gas_temperature.value))]


aul = nacltbl['Aij'].quantity
nu_ul = nacltbl['Frequency'].quantity
eu = nacltbl['E_U'].quantity
el = nacltbl['E_L'].quantity
gu = nacltbl['gu'].quantity
gl = nacltbl['gl'].quantity


oscillator_strength_ul = (aul / (2*np.pi*(nu_ul)**2*constants.e.si**2 /
                                 (constants.eps0 * constants.m_e *
                                  constants.c**3) * (gu/gl)).decompose())
bul = (constants.e.si**2 / (4*constants.eps0*constants.m_e*constants.h*nu_ul) *
       oscillator_strength_ul).decompose()
blu = bul * gl/gu

if 'Bij' not in nacltbl.colnames:
    nacltbl.add_column(Column(name='Bij', data=bul))
    nacltbl.add_column(Column(name='Bji', data=blu))
    nacltbl.add_column(Column(name='Fij', data=oscillator_strength_ul))

# 4 pi / c * B_nu
radiation_field_energy_density = (4*np.pi/constants.c *
                                  (2*constants.h*nu_ul**3/constants.c**2) *
                                  (np.exp(constants.h*nu_ul /
                                          (constants.k_B * stellar_temperature)) -
                                   1)**-1)
# Dilute the radiation field by 1/r^2
radiation_field_energy_density /= (35*u.au / stellar_radius)**2


# Compute the collisional excitation, which will be applied only for v=0
nu = (gu * np.exp(-eu/gas_temperature) / Zgas).decompose()
nl = (gl * np.exp(-el/gas_temperature) / Zgas).decompose()

#nl = ((aul * nu + bul * nu * radiation_field_energy_density) /
#      (blu * radiation_field_energy_density)).decompose()
nl_nu_eq = nl/nu
nl_nu_rad = ((aul + bul * radiation_field_energy_density) /
             (blu * radiation_field_energy_density)).decompose()

# energy radiated divided by 4pi steradians
# per Hz because... defined that way?  Will need to integrate over line profile...
# This integrated over a path gives the specific intensity (Jy)
emission_coeff = constants.h * nu_ul / (4*np.pi*u.sr) * nu * voldens * (aul + bul * radiation_field_energy_density) / u.Hz
if 'emi_coeff_4000K' not in nacltbl.colnames:
    nacltbl.add_column(Column(name='emi_coeff_4000K', data=emission_coeff))


sel = (#(nacltbl['vu'] == 2) &
       #(nacltbl['vl'] == 1) &
       (nacltbl['vu'] < 20) &
       (nacltbl['vu'] == nacltbl['vl']) &
       (nacltbl['Ju'] == 18) &
       (nacltbl['Jl'] == 17))

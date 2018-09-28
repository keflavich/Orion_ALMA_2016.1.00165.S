import numpy as np
import paths
import parse_exomol_files
from astropy import units as u
from astropy import constants
from astropy.table import Table

if 'nacl_exomol' not in locals():
    nacl_exomol = parse_exomol_files.ExoMol('NaCl', '23Na-35Cl')

nacltbl = nacl_exomol.transitions


nacltbl[(nacltbl['vu'] == 2) & (nacltbl['vl'] == 1) & (nacltbl['Ju'] == 18)]
nacltbl[(nacltbl['vu'] == 2) & (nacltbl['vl'] == 2) & (nacltbl['Ju'] == 18)]

sel = ((nacltbl['vu'] == 2) &
       (nacltbl['vl'] == 1) &
       (nacltbl['Ju'] == 18) &
       (nacltbl['Jl'] == 17))
aul = nacltbl[sel]['Aij'].quantity
nuul = nacltbl[sel]['Frequency'].quantity
eu = nacltbl[sel]['E_U'].quantity
el = nacltbl[sel]['E_L'].quantity
gu = nacltbl[sel]['gu'].quantity
gl = nacltbl[sel]['gl'].quantity


oscillator_strength_ul = (aul / (2*np.pi*(nuul)**2*constants.e.si**2 /
                                 (constants.eps0 * constants.m_e *
                                  constants.c**3) * (gu/gl)).decompose())
bul = (constants.e.si**2 / (4*constants.eps0*constants.m_e*constants.h*nuul) *
       oscillator_strength_ul).decompose()
blu = bul * gl/gu

print(aul, bul, blu, oscillator_strength_ul)

stellar_temperature = 4000*u.K
stellar_luminosity = 2e4*u.L_sun
stellar_radius = ((stellar_luminosity / (constants.sigma_sb *
                                         stellar_temperature**4) /
                   (4*np.pi))**0.5).to(u.R_sun)
# 4 pi / c * B_nu
radiation_field_energy_density = (4*np.pi/constants.c *
                                  (2*constants.h*nuul**3/constants.c**2) *
                                  (np.exp(constants.h*nuul /
                                          (constants.k_B * stellar_temperature)) -
                                   1)**-1)
# Dilute the radiation field by 1/r^2
radiation_field_energy_density /= (35*u.au / stellar_radius)**2

# partition function: look up from table
pf = Table.read(paths.salty('23Na-37Cl__Barton.pf'), format='ascii')

Z = pf['col2'][np.argmin(np.abs(pf['col1'] - stellar_temperature.value))]

nu = (gu * np.exp(-eu/stellar_temperature) / Z).decompose()
nl = (gl * np.exp(-el/stellar_temperature) / Z).decompose()
#nl = ((aul * nu + bul * nu * radiation_field_energy_density) /
#      (blu * radiation_field_energy_density)).decompose()
nl_nu_eq = nl/nu
nl_nu_rad = ((aul + bul * radiation_field_energy_density) /
             (blu * radiation_field_energy_density)).decompose()

print(nl_nu_rad, nl_nu_eq)

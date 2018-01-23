"""
Determine the offset velocity for each of the lines
"""
import numpy as np
import pyspeckit
import lines
import paths

from astroquery.splatalogue import Splatalogue

from astropy import table
from astropy import units as u
from astropy import constants


dv = 15 * u.km/u.s
v = 5.5 * u.km/u.s
dv_linesearch = 2.5*u.km/u.s

fits = {}

for spw in (0,1,2,3):
    sp = pyspeckit.Spectrum(paths.dpath('stacked_spectra/OrionSourceI_B6_spw{0}.fits'
                                        .format(spw)),)

    for linename, freq in lines.disk_lines.items():

        xmin = freq*(1+(v-dv)/constants.c)
        xmax = freq*(1+(v+dv)/constants.c)

        slc = sp.slice(xmin,xmax)
        if len(slc) == 0:
            continue


        guesses = [slc.data.max(),
                   (freq*(1+v/constants.c)).to(u.GHz).value,
                   (2*u.km/u.s/constants.c*freq).to(u.GHz).value]
        #print(guesses)

        slc.specfit(guesses=guesses)
        #print(slc.specfit.parinfo)

        frq = u.Quantity(slc.specfit.parinfo['SHIFT0'], u.GHz)
        result = Splatalogue.query_lines(frq + (v-dv_linesearch)/constants.c*frq,
                                         frq + (v+dv_linesearch)/constants.c*frq)
        ref = np.array(result['Freq-GHz'])*u.GHz
        result.add_column(table.Column(name='velocity', data=-((frq-ref)/(ref) * constants.c).to(u.km/u.s)))
        linesearch = result['Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','velocity', 'E_U (K)']

        fits[linename] = {'pars': slc.specfit.parinfo,
                          'vel':
                          ((u.Quantity(slc.specfit.parinfo['SHIFT0'].value,
                                       u.GHz) - freq) / freq *
                           constants.c.to(u.km/u.s)),
                          'vwidth':
                          ((u.Quantity(slc.specfit.parinfo['WIDTH0'].value,
                                       u.GHz)) / freq *
                           constants.c.to(u.km/u.s)),
                          'linesearch': linesearch
                         }


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

import latex_info

dv = 15 * u.km/u.s
v = 5.5 * u.km/u.s
dv_linesearch = 2.5*u.km/u.s

fits = {}

bad_fits = ['Unknown_8', # only half the line is detected
           ]

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
                          'linesearch': linesearch,
                          'freq': freq,
                          'spectrum': slc,
                         }

linenames = table.Column(name='Line Name', data=sorted(fits.keys()))
freqs = table.Column(name='Frequency', data=u.Quantity([fits[ln]['freq'] for ln in linenames]))
velos = table.Column(name='Fitted velocity', data=u.Quantity([fits[ln]['vel'] for ln in linenames]))
vwidths = table.Column(name='Fitted Width', data=u.Quantity([fits[ln]['vwidth'] for ln in linenames]))
ampls = table.Column(name='Fitted Amplitude', data=u.Quantity([fits[ln]['pars']['AMPLITUDE0'].value*1e3 for ln in linenames], u.mJy))

tbl = table.Table([linenames, freqs, velos, vwidths, ampls])

tbl.write('fitted_stacked_lines.txt', format='ascii.fixed_width')

badmask = np.array([ln in bad_fits for ln in linenames], dtype='bool')

linenames = table.Column(["U{0:0.3f}".format(freq) if 'Unknown' in ln else ln
                          for ln, freq in zip(linenames, freqs)],
                         name='Line Name',
                        )


tbl = table.Table([linenames, freqs, velos, vwidths, ampls])
tbl['Fitted Width'][badmask] = np.nan
tbl['Fitted Amplitude'][badmask] = np.nan
ulines = np.array(['U' in ln for ln in linenames], dtype='bool')
tbl = tbl[ulines]['Line Name', 'Frequency', 'Fitted Width', 'Fitted Amplitude']

tbl.sort('Frequency')

tbl.write('unknown_line_fits.txt', format='ascii.fixed_width')


latexdict = latex_info.latexdict.copy()
latexdict['header_start'] = '\label{tab:unknown_line_frequencies}'
latexdict['caption'] = 'Unknown Line Frequencies'
latexdict['preamble'] = '\centering'
latexdict['tablefoot'] = ('\n\par The frequencies listed have a systematic '
                          'uncertainty of about 2 \\kms (1.5 MHz) because '
                          'they are referenced to the U232.511 line, which '
                          'has an unknown rest frequency.  The rest frequency '
                          'used for the U232.511 line was selected to maximize '
                          'the symmetry of the emission around 5 \\kms.  '
                          'Some lines were detected in only part of the disk '
                          'and therefore had bad or malformed profiles in the '
                          'stacked spectrum; these have fits marked with -\'s.'
                         )
formats = {'Frequency': lambda x: "{0:0.3f}".format(x),
           'Fitted Width': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
          }
tbl.write(paths.texpath('unknown_line_freqs.tex'),
          formats=formats,
          latexdict=latexdict,
          overwrite=True)

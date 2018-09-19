"""
Determine the line parameters for each of the lines
"""
import numpy as np
import pyspeckit
import lines
import paths

from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
from astropy.io import fits

from astropy import table
from astropy import stats
from astropy import units as u
from astropy import constants

import pylab as pl

from salt_tables import salt_tables

import latex_info

dv = 15 * u.km/u.s
v = 5.5 * u.km/u.s
dv_linesearch = 10.0*u.km/u.s

linefits = {}

bad_fits = ['Unknown_8', # only half the line is detected
           ]

chem_re = "KCl|NaCl|K37Cl|Na37Cl"


if 'doplot' not in locals():
    doplot = False

if doplot:
    pl.figure(0).clf()

for spw in (0,1,2,3):
    for band in ('B3', 'B6', 'B7.lb'):
        fn = paths.dpath('stacked_spectra/OrionSourceI_{band}_spw{0}_robust0.5.fits'
                         .format(spw, band=band))
        sp = pyspeckit.Spectrum(fn)

        rms = stats.mad_std(sp.data)
        sp.error[:] = rms
        print(rms)

        beams = fits.open(fn)[1]
        beam_area = np.median(beams.data['BMAJ'] * beams.data['BMIN'] * np.pi *
                              u.arcsec**2)
        jtok = u.brightness_temperature(frequency=sp.xarr.mean(),
                                        beam_area=beam_area)

        for linename, freq in lines.disk_lines.items():

            xmin = freq*(1+(v-dv)/constants.c)
            xmax = freq*(1+(v+dv)/constants.c)

            slc = sp.slice(xmin,xmax)
            if len(slc) == 0:
                continue


            guesses = [np.max([slc.data.max(), 0.05]),
                       (freq*(1+v/constants.c)).to(u.GHz).value,
                       (2*u.km/u.s/constants.c*freq).to(u.GHz).value]
            #print(guesses)

            if doplot:
                slc.plotter(figure=pl.figure(0), clear=True)

            slc.specfit(guesses=guesses,
                        limits=[(0,1),
                                (xmin.value, xmax.value),
                                (0, 15)],
                        limited=[(True,True)]*3,
                       )

            if doplot:
                slc.plotter.savefig(paths.fpath('spectral_fits/{linename}_{band}_{spw}_{freq}.png'
                                                .format(linename=linename,
                                                        band=band,
                                                        freq=freq,
                                                        spw=spw)))

            frq = u.Quantity(slc.specfit.parinfo['SHIFT0'], u.GHz)
            result = Splatalogue.query_lines(freq - (dv_linesearch)/constants.c*freq,
                                             freq + (dv_linesearch)/constants.c*freq,
                                             chemical_name=chem_re
                                            )
            for tbl in salt_tables:
                match = ((tbl['Freq'].quantity > freq - (dv_linesearch)/constants.c*freq) &
                         (tbl['Freq'].quantity < freq + (dv_linesearch)/constants.c*freq))
                result = tbl[match]
                #print(match.any())
                if match.any():
                    #print("Matched {0}".format(linename))
                    break

            #if len(result) > 0:
            #    result = mt(result)
            if len(result) >= 1:
                if len(result) > 1:
                    print(result)
                #eu = u.Quantity(result[0]['E_U (cm^-1)'], u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()).value
                eu = u.Quantity(result[0]['E_U'], u.K)
                species = result[0]['Species']
                #qn = result[0]['Resolved QNs']
                qn = result[0]['QNs']
                aul = result[0]['Aij']
                deg = result[0]['gu']
            else:
                #print("No match for {0}".format(linename))
                eu = u.Quantity(np.nan, u.K)
                species = ''
                qn = ''
                aul = 0
                deg = 0

            #ref = np.array(result['Freq'])*u.GHz
            #result.add_column(table.Column(name='velocity', data=-((frq-ref)/(ref) * constants.c).to(u.km/u.s)))
            linesearch = result#['Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','velocity', 'E_U (K)']

            linefits[linename] = {'pars': slc.specfit.parinfo,
                                  'vel':
                                  ((u.Quantity(slc.specfit.parinfo['SHIFT0'].value,
                                               u.GHz) - freq) / freq *
                                   constants.c.to(u.km/u.s)),
                                  'evel':
                                  ((u.Quantity(slc.specfit.parinfo['SHIFT0'].error,
                                               u.GHz)) / freq *
                                   constants.c.to(u.km/u.s)),
                                  'vwidth':
                                  ((u.Quantity(slc.specfit.parinfo['WIDTH0'].value,
                                               u.GHz)) / freq *
                                   constants.c.to(u.km/u.s)),
                                  'evwidth':
                                  ((u.Quantity(slc.specfit.parinfo['WIDTH0'].error,
                                               u.GHz)) / freq *
                                   constants.c.to(u.km/u.s)),
                                  'linesearch': linesearch,
                                  'freq': freq,
                                  'spectrum': slc,
                                  'EU_K': eu,
                                  'species': species,
                                  'qn': qn,
                                  'jtok': (1*u.Jy).to(u.K, equivalencies=jtok),
                                  'aul': aul,
                                  'deg': deg,
                                 }

linenames = table.Column(name='Line Name', data=sorted(linefits.keys()))
freqs = table.Column(name='Frequency', data=u.Quantity([linefits[ln]['freq'] for ln in linenames]))
velos = table.Column(name='Fitted Velocity', data=u.Quantity([linefits[ln]['vel'] for ln in linenames]))
vwidths = table.Column(name='Fitted Width', data=u.Quantity([linefits[ln]['vwidth'] for ln in linenames]))
evelos = table.Column(name='Fitted Velocity error', data=u.Quantity([linefits[ln]['evel'] for ln in linenames]))
evwidths = table.Column(name='Fitted Width error', data=u.Quantity([linefits[ln]['evwidth'] for ln in linenames]))
ampls = table.Column(name='Fitted Amplitude', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].value*1e3 for ln in linenames], u.mJy))
amplsK = table.Column(name='Fitted Amplitude K', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].value*linefits[ln]['jtok'].value for ln in linenames], u.K))
eampls = table.Column(name='Fitted Amplitude error', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].error*1e3 for ln in linenames], u.mJy))
eamplsK = table.Column(name='Fitted Amplitude error K', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].error*linefits[ln]['jtok'].value for ln in linenames], u.K))
jtok = table.Column(name='Jy/K', data=u.Quantity([linefits[ln]['jtok'].value for ln in linenames], u.Jy/u.K))
eu = table.Column(name='EU_K', data=u.Quantity([linefits[ln]['EU_K'] for ln in linenames], u.K))
species = table.Column(name='Species', data=[linefits[ln]['species'] for ln in linenames])
qn = table.Column(name='QNs', data=[linefits[ln]['qn'] for ln in linenames])
deg = table.Column(name='deg', data=[linefits[ln]['deg'] for ln in linenames])
Aij = table.Column(name='Aij', data=[linefits[ln]['aul'] for ln in linenames])

tbl1 = table.Table([linenames, species, qn, freqs, velos, evelos, vwidths, evwidths, ampls, eampls, amplsK, eamplsK, jtok, eu, deg, Aij])

tbl1.write(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

badmask = np.array([ln in bad_fits for ln in linenames], dtype='bool')

linenames = table.Column(["U{0:0.3f}".format(freq) if 'Unknown' in ln else ln
                          for ln, freq in zip(linenames, freqs)],
                         name='Line Name',
                        )


tbl = table.Table([linenames, species, qn, freqs, velos, vwidths, ampls, eu])
tbl['Fitted Width'][badmask] = np.nan
tbl['Fitted Amplitude'][badmask] = np.nan
ulines = np.array(['U' in ln for ln in linenames], dtype='bool')
tbl = tbl[ulines]['Line Name', 'Species', 'QNs', 'Frequency', 'Fitted Width', 'Fitted Amplitude', 'EU_K']

tbl.sort('Frequency')

tbl.write(paths.tpath('line_fits.txt'), format='ascii.fixed_width')


latexdict = latex_info.latexdict.copy()
latexdict['header_start'] = '\label{tab:unknown_line_frequencies}'
latexdict['caption'] = 'Unknown Line Frequencies'
latexdict['preamble'] = '\centering'
latexdict['tablefoot'] = ('\n\par '
                         )
formats = {'Frequency': lambda x: "{0:0.3f}".format(x),
           'Fitted Width': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
          }
tbl.write(paths.texpath2('line_parameters.tex'),
          formats=formats,
          latexdict=latexdict,
          overwrite=True)



# plot some things....

pl.figure(1).clf()

kclmask = np.array(['KCl' in row['Species'] for row in tbl1])
k41clmask = np.array(['41KCl' in row['Species'] for row in tbl1])
k37clmask = np.array(['K37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][kclmask], tbl1['Fitted Amplitude'][kclmask], 'o', label='KCl')
pl.plot(tbl1['EU_K'][k37clmask], tbl1['Fitted Amplitude'][k37clmask], 's', label='K$^{37}$Cl')
pl.plot(tbl1['EU_K'][k41clmask], tbl1['Fitted Amplitude'][k41clmask], 'd', label='$^{41}$KCl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [mJy]")
pl.legend(loc='best')
pl.savefig(paths.fpath('KCl_amp_vs_eu.pdf'))

pl.figure(2).clf()
naclmask = np.array(['NaCl' in row['Species'] for row in tbl1])
na37clmask = np.array(['Na37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][naclmask], tbl1['Fitted Amplitude'][naclmask], 'o', label='NaCl')
pl.plot(tbl1['EU_K'][na37clmask], tbl1['Fitted Amplitude'][na37clmask], 's', label='Na$^{37}$Cl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [mJy]")
pl.legend(loc='best')
pl.savefig(paths.fpath('NaCl_amp_vs_eu.pdf'))

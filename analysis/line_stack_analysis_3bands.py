"""
Determine the line parameters for each of the lines
"""
import re
import numpy as np
import pyspeckit
import lines
import paths

from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
from astropy.io import fits
from astropy.table import Table, Column

from astropy import table
from astropy import stats
from astropy import units as u
from astropy import constants

import pylab as pl

from salt_tables import KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37

import latex_info

# Filter out the high-v lines; they result in confusion in some cases (there
# are some v=10 lines really close to v=8 lines...)
salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]
salt_tables = [tbl[tbl['vu']<9] for tbl in salt_tables]

dv = 15 * u.km/u.s
v = 5.5 * u.km/u.s
dv_linesearch = 5.0*u.km/u.s

linefits = {}

chem_re = "KCl|NaCl|K37Cl|Na37Cl"

detection_table = Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')
nondetections = (detection_table['Flag'] == '-n') | (detection_table['Flag'] == 'cn')
detection_table = detection_table[~nondetections]

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


        #for linename, freq in lines.disk_lines.items():
        for row in detection_table:

            linename = row['Species']
            freq = u.Quantity(row['Frequency'], u.GHz)
            detection = row['Flag'][1] == 'd'
            if not detection:
                continue

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
                                  'flag': row['Flag'],
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
flag = table.Column(name='Flag', data=[linefits[ln]['flag'] for ln in linenames])

vre = re.compile('v=([0-9]+)')
vstate = [int(vre.search(ss).groups()[0]) for ss in species]
jre = re.compile('J=([0-9]+)-([0-9]+)')
Ju,Jl = zip(*[map(int, jre.search(ss).groups()) for ss in species])
qnv = (Column(name='v', data=vstate))
qnju = (Column(name='J$_u$', data=Ju))
qnjl = (Column(name='J$_l$', data=Jl))


tbl1 = table.Table([linenames, species, qn, qnv, qnju, qnjl, freqs, velos, evelos, vwidths, evwidths, ampls, eampls, amplsK, eamplsK, jtok, eu, deg, Aij, flag, ])

tbl1.write(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')



# create a subtable for the paper

linenames = table.Column([lines.texnames[ln]
                          if ln in lines.texnames
                          else ln
                          for ln, freq in zip(linenames, freqs)
                         ],
                         name='Line Name',
                        )


tbl = table.Table([linenames, qnv, qnju, qnjl, freqs, velos, evelos, vwidths, evwidths, amplsK, eamplsK, eu, flag])

tbl.sort('Frequency')

bad_fits = []

badmask = np.array([ln in bad_fits for ln in linenames], dtype='bool')
badmask |= ((tbl['Fitted Width error'] > tbl['Fitted Width']) |
            (tbl['Fitted Velocity error'] > 5) |
            np.array([flg[1] in 'nq' for flg in tbl['Flag']])
           )


tbl.write(paths.tpath('line_fits.txt'), format='ascii.fixed_width')

formats = {'Frequency': lambda x: "{0:0.5f}".format(x),
           'Fitted Width': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Width error': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Velocity': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Velocity error': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude error K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'EU_K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
          }
rename = {'Fitted Width':'Width',
          'Fitted Width error':'Width error',
          'Fitted Velocity':'Velocity',
          'Fitted Velocity error':'Velocity error',
          'Fitted Amplitude error K':'Amplitude error',
          'Fitted Amplitude K':'Amplitude',
          'EU_K': 'E$_U$',
         }
maskcols = ['Fitted Width',
            'Fitted Width error',
            'Fitted Velocity',
            'Fitted Velocity error',
            'Fitted Amplitude error K',
            'Fitted Amplitude K', ]

for msk in maskcols:
    tbl[msk][badmask] = np.nan
    tbl[badmask][msk] = np.nan

for old, new in rename.items():
    if old in tbl.columns:
        tbl.rename_column(old, new)
    formats[new] = formats[old]

#print(tbl)


def label_by_vstate(tablepath, ncols):
    with open(tablepath, 'r') as fh:
        lines = fh.readlines()

    started_data = False
    current_vstate = None
    vstate_lines = []
    with open(tablepath, 'w') as fh:
        for line in lines:
            if line.startswith(r'\hline'):
                started_data = True
            elif line.startswith(r'\end{tabular}'):
                started_data = False

                headerstring = (r"&\vspace{{-0.75em}}\\""\n"
                                r"\multicolumn{{{ncol}}}{{c}}{{$v = {vstate}$}} \\""\n"
                                r"\vspace{{-0.75em}}\\""\n").format(ncol=ncols-1,
                                                                    vstate=current_vstate)
                if len(vstate_lines) > 0:
                    fh.write(headerstring)
                    for ll in vstate_lines:
                        fh.write(ll)
            elif line.startswith('v &') or (not started_data and line.startswith(r' &')):
                # this is the header line or unit line
                fh.write("&".join(line.split("&")[1:]))
                continue
            elif line.startswith(r'\begin{tabular}'):
                # strip out one column
                # (last character is \n)
                line = line[:-3] + line[-2:]
                fh.write(line)
                continue
            elif line[0].isdigit() and started_data:
                vstate = int(line[0])
                if current_vstate is None:
                    current_vstate = vstate
                    line = "&".join(line.split("&")[1:])
                    vstate_lines.append(line)
                elif current_vstate == vstate:
                    # drop v= part
                    line = "&".join(line.split("&")[1:])
                    vstate_lines.append(line)
                else:
                    line = "&".join(line.split("&")[1:])
                    vstate_lines.append(line)
                    headerstring = (r"&\vspace{{-0.75em}}\\""\n"
                                    r"\multicolumn{{{ncol}}}{{c}}{{$v = {vstate}$}} \\""\n"
                                    r"\vspace{{-0.75em}}\\""\n").format(ncol=ncols-1,
                                                                        vstate=current_vstate)
                    if len(vstate_lines) > 0:
                        fh.write(headerstring)
                        for ll in vstate_lines:
                            fh.write(ll)
                    else:
                        print("For v={0}, no lines are found.".format(current_vstate))

                    # reset
                    vstate_lines = []
                    current_vstate = vstate
                continue

            fh.write(line)

def merge_errors(tbl, datacol, errorcol):
    new_col = ["{0} ({1})".format(formats[datacol](row[datacol]),
                                  formats[errorcol](row[errorcol]))
               for row in tbl]
    tbl.rename_column(datacol, '_'+datacol)
    tbl.add_column(Column(data=new_col, name=datacol, unit=tbl['_'+datacol].unit))

merge_errors(tbl, 'Width', 'Width error')
merge_errors(tbl, 'Amplitude', 'Amplitude error')
merge_errors(tbl, 'Velocity', 'Velocity error')
del formats['Velocity']
del formats['Width']
del formats['Amplitude']

tbl = tbl['Line Name', 'v', 'J$_u$', 'J$_l$', 'Frequency', 'Velocity',
          'Width', 'Amplitude', 'E$_U$', ]


#for salt in ('NaCl', 'Na$^{37}Cl', 'KCl', 'K$^{37}$Cl', '$^{41}$KCl',
#             '$^{41}$K$^{37}$Cl'):
salt_to_barton = {'NaCl': '23Na-35Cl',
                  'Na37Cl': '23Na-37Cl',
                  'KCl': '39K-35Cl',
                  'K37Cl': '39K-37Cl',
                  '41KCl': '41K-35Cl',
                  '41K37Cl': '41K-37Cl'}
for salt in ('NaCl', 'Na37Cl', 'KCl', 'K37Cl', '41KCl', '41K37Cl'):
    texsalt = salt.replace("37","$^{37}$").replace("41","$^{41}$")
    latexdict = latex_info.latexdict.copy()
    latexdict['header_start'] = '\label{{tab:{0}_salt_lines}}'.format(salt)
    latexdict['caption'] = '{0} Lines'.format(texsalt)
    latexdict['preamble'] = '\centering'
    latexdict['tablefoot'] = ('\n\par '
                             )
    tbl.sort('v')
    mask = np.array([ln.startswith(salt_to_barton[salt]) for ln in tbl['Line Name']])
    mask &= ~badmask
    print("{0} matches for {1}".format(mask.sum(), salt))
    for row in tbl[mask]:
        assert row['Line Name'].startswith(salt_to_barton[salt])
    columns = tbl.colnames[1:] # drop Line Name
    tbl[mask][columns].write(paths.texpath2('{0}_line_parameters.tex'.format(salt)),
                             formats=formats,
                             latexdict=latexdict,
                             overwrite=True)
    #print(tbl[mask][columns])
    assert len(tbl[mask][columns]) == mask.sum()
    label_by_vstate(paths.texpath2('{0}_line_parameters.tex'.format(salt)),
                    ncols=len(columns)
                   )



# plot some things....

pl.figure(1).clf()

kclmask = np.array(['K-35Cl' in row['Species'] for row in tbl1])
k41clmask = np.array(['41K-35Cl' in row['Species'] for row in tbl1])
k37clmask = np.array(['K-37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][kclmask], tbl1['Fitted Amplitude'][kclmask], 'o', label='KCl')
pl.plot(tbl1['EU_K'][k37clmask], tbl1['Fitted Amplitude'][k37clmask], 's', label='K$^{37}$Cl')
pl.plot(tbl1['EU_K'][k41clmask], tbl1['Fitted Amplitude'][k41clmask], 'd', label='$^{41}$KCl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [mJy]")
pl.legend(loc='best')
pl.savefig(paths.fpath('KCl_amp_vs_eu.pdf'))

pl.figure(2).clf()
naclmask = np.array(['Na-35Cl' in row['Species'] for row in tbl1])
na37clmask = np.array(['Na-37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][naclmask], tbl1['Fitted Amplitude'][naclmask], 'o', label='NaCl')
pl.plot(tbl1['EU_K'][na37clmask], tbl1['Fitted Amplitude'][na37clmask], 's', label='Na$^{37}$Cl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [mJy]")
pl.legend(loc='best')
pl.savefig(paths.fpath('NaCl_amp_vs_eu.pdf'))

"""
Determine the line parameters for each of the lines
"""
import re
import numpy as np
import pyspeckit
import lines
import paths
import glob

from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
from astropy.io import fits
from astropy.table import Table, Column
from astropy import coordinates

from astropy import table
from astropy import stats
from astropy import units as u
from astropy import constants

import pylab as pl
from spectral_cube import SpectralCube

from salt_tables import KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37

import latex_info

# Filter out the high-v lines; they result in confusion in some cases (there
# are some v=10 lines really close to v=8 lines...)
salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]
salt_tables = [tbl[tbl['vu']<9] for tbl in salt_tables]

dv = 20 * u.km/u.s
vcen = v = 5.5 * u.km/u.s
dv_linesearch = 5.0*u.km/u.s

topright = coordinates.SkyCoord('5:35:14.5134090216', '-5:22:30.5431496522', frame='icrs', unit=(u.h, u.deg))
trv = 19.0*u.km/u.s
bottomleft = coordinates.SkyCoord('5:35:14.5211504272', '-5:22:30.7054489346', frame='icrs', unit=(u.h, u.deg))
blv = -7.6*u.km/u.s

linefits = {}

chem_re = "KCl|NaCl|K37Cl|Na37Cl"

detection_table = Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')
nondetections = (detection_table['Flag'] == '-n') | (detection_table['Flag'] == 'cn')
detection_table = detection_table[~nondetections]

if 'doplot' not in locals():
    doplot = False

if doplot:
    pl.figure(0).clf()

epath = '/orange/adamginsburg/orion/2016.1.00165.S/external/'

spwband = {spw: {band: f'{epath}/OrionSourceI_only.{band}.robust0.5.spw{spw}.*clarkclean10000_medsub.image.pbcor.fits'
    for band in ('B3', 'B6', 'B7.lb')}
    for spw in (0,1,2,3)}
for spw in (0,1,2):
    spwband[spw]['B8_495'] = f'/orange/adamginsburg/salt/hirota/ALMA-SourceI/band8_495G_spw{spw}-subim.FITS'
for spw in (0,1,3):
    spwband[spw]['B10_850'] = f'/orange/adamginsburg/salt/hirota/ALMA-SourceI/band10_850G_spw{spw}-subim.FITS'

for spw in (0,1,2,3):
    for band in ('B10_850', 'B3', 'B6', 'B7.lb', 'B8_495', ):
        if band in spwband[spw]:
            fn = spwband[spw][band]
        else:
            continue
        flist = glob.glob(fn)
        if len(flist) == 0:
            print(f"Skipped spw {spw} band {band}")
            continue
        assert len(flist) > 0
        cube = SpectralCube.read(flist[0]).with_spectral_unit(u.GHz)
        if hasattr(cube, 'beams'):
            jtok = cube.beams.common_beam().jtok(cube.spectral_axis).mean()
        else:
            jtok = cube.beam.jtok(cube.spectral_axis).mean()

        for cen, vel, cname in ((topright, trv, 'topright'), (bottomleft, blv, 'bottomleft')): 

            xx,yy = map(int, cube.wcs.celestial.world_to_pixel(cen))
            csp = cube[:,yy,xx].to(u.K)

            sp = pyspeckit.Spectrum(data=csp.value, unit=csp.unit, xarr=csp.spectral_axis)
            sp.xarr.convert_to_unit(u.GHz)

            rms = stats.mad_std(sp.data)
            sp.error[:] = rms
            print(f"RMS={rms}")

            if doplot:
                sp.plotter(figure=pl.figure(1, figsize=(16,6)), clear=True)

            linenames, linefreqs = [],[]
            #for linename, freq in lines.disk_lines.items():
            for row in detection_table:

                linename = row['Species']
                freq = u.Quantity(row['Frequency'], u.GHz)
                detection = (row['Flag'][1] == 'd') or (row['Flag'] == '--')
                if not detection:
                    continue

                xmin = (vel - dv).to(u.GHz, u.doppler_radio(freq))
                xmax = (vel + dv).to(u.GHz, u.doppler_radio(freq))

                slc = sp.slice(xmin,xmax)
                if len(slc) == 0:
                    continue

                slc.xarr.convert_to_unit(u.km/u.s, refX=freq, velocity_convention='radio')
                print(f"{cname} vel={vel} xmin/xmax: {xmin},{xmax}, {slc.xarr.min()}, {slc.xarr.max()}")

                guesses = [np.max([slc.data.max(), 0.05]),
                           (freq*(1+vel/constants.c)).to(u.GHz).value,
                           (2*u.km/u.s/constants.c*freq).to(u.GHz).value]
                guesses = [np.max([slc.data.max(), 25]),
                           vel.value,
                           2]
                print(f"{linename}: guesses={guesses}")

                if doplot:
                    slc.plotter(figure=pl.figure(0), clear=True)

                slc.specfit(guesses=guesses,
                            limits=[(0,500),
                                    (-40,50),#(xmin.value, xmax.value),
                                    (0, 10)],
                            limited=[(True,True)]*3,
                            verbose=False,
                           )

                if doplot:
                    sp.plotter.axis.plot(slc.xarr.as_unit(u.GHz), slc.specfit.model, color='r', linewidth=0.5)
                    slc.plotter.savefig(paths.fpath(f'spectral_fits/{linename}_{band}_{spw}_{freq.value}_{cname}.png'))

                fitvel = u.Quantity(slc.specfit.parinfo['SHIFT0'], u.km/u.s)
                frq = fitvel.to(u.GHz, u.doppler_radio(freq))
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

                vfit = fitvel
                if slc.specfit.parinfo['SHIFT0'].error is None:
                    evel = np.nan*u.km/u.s
                else:
                    evel = u.Quantity(slc.specfit.parinfo['SHIFT0'].error, u.km/u.s)
                vwidth = u.Quantity(slc.specfit.parinfo['WIDTH0'].value, u.km/u.s)
                if slc.specfit.parinfo['WIDTH0'].error is None:
                    evwidth = np.nan*u.km/u.s
                else:
                    evwidth = ((u.Quantity(slc.specfit.parinfo['WIDTH0'].error,
                                                       u.GHz)) / freq *
                                           constants.c.to(u.km/u.s))
                try:
                    SNR = slc.specfit.parinfo['AMPLITUDE0'].value / slc.specfit.parinfo['AMPLITUDE0'].error
                except TypeError:
                    SNR = 0
                
                if SNR > 0:
                    linenames.append(linename)
                    linefreqs.append(freq)
                    print(f"Line {linename} with frequency {freq} has SNR={SNR:0.1f}")

                # initialization
                if linename not in linefits:
                    linefits[linename] = {}

                linefits[linename][cname] = {'pars': slc.specfit.parinfo,
                                      'vel': vfit,
                                      'evel': evel,
                                      'vwidth': vwidth,
                                      'evwidth': evwidth,
                                      'linesearch': linesearch,
                                      'freq': freq,
                                      'spectrum': slc,
                                      'EU_K': eu,
                                      'species': species,
                                      'qn': qn,
                                      'jtok': jtok,
                                      'aul': aul,
                                      'deg': deg,
                                      'snr': SNR,
                                      'flag': row['Flag'],
                                      'corner': cname,
                                     }
            if doplot:
                sp.plotter.axis.set_ylim(-20, 250)
                sp.plotter.line_ids(linenames, linefreqs, velocity_offset=vel,
                                       label1_size=16,
                                       auto_yloc_fraction=0.75)
                for txt in sp.plotter.axis.texts:
                    txt.set_backgroundcolor((1,1,1,0.9))
                for obj in sp.plotter.axis.texts+sp.plotter.axis.lines:
                    if 'Na' in obj.get_label():
                        obj.set_color('r')
                        obj.set_zorder(5)
                    elif 'K' in obj.get_label():
                        obj.set_color('b')
                        obj.set_zorder(10)
                sp.plotter.axis.set_ylim(-20, 250)

                sp.plotter.savefig(paths.fpath(f'spectral_fits/SrcI_{band}_{spw}_{cname}_fits.png'))

linenames = table.Column(name='Line Name', data=sorted(linefits.keys()))
def makecol(colname, unit=None, linefits=linefits, parname=None, error=False, linenames=linenames):
    if parname is not None:
        if error:
            return u.Quantity([linefits[ln][cname][colname][parname].error
                if linefits[ln][cname][colname][parname].error else np.nan
                for ln in linenames for cname in linefits[ln]], unit)
        else:
            return u.Quantity([linefits[ln][cname][colname][parname].value for ln in linenames for cname in linefits[ln]], unit)
    if unit is not None:
        return u.Quantity([linefits[ln][cname][colname] for ln in linenames for cname in linefits[ln]], unit)
    else:
        return [linefits[ln][cname][colname] for ln in linenames for cname in linefits[ln]]

freqs = table.Column(name='Frequency', data=makecol('freq', unit=u.GHz))
velos = table.Column(name='Fitted Velocity', data=makecol('vel', u.km/u.s))
vwidths = table.Column(name='Fitted Width', data=makecol('vwidth', u.km/u.s))
evelos = table.Column(name='Fitted Velocity error', data=makecol('evel', unit=u.km/u.s))
evwidths = table.Column(name='Fitted Width error', data=makecol('evwidth', unit=u.km/u.s))
#ampls = table.Column(name='Fitted Amplitude', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].value*1e3 for ln in linenames], u.mJy))
#amplsK = table.Column(name='Fitted Amplitude K', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].value*linefits[ln]['jtok'].value for ln in linenames], u.K))
amplsK = table.Column(name='Fitted Amplitude K', data=makecol('pars', parname='AMPLITUDE0', unit=u.K))
#eampls = table.Column(name='Fitted Amplitude error', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].error*1e3 for ln in linenames], u.mJy))
#eamplsK = table.Column(name='Fitted Amplitude error K', data=u.Quantity([linefits[ln]['pars']['AMPLITUDE0'].error*linefits[ln]['jtok'].value for ln in linenames], u.K))
eamplsK = table.Column(name='Fitted Amplitude error K', data=makecol('pars', parname='AMPLITUDE0', unit=u.K, error=True))
integrated = table.Column(name='Integrated Intensity', data=amplsK.quantity*vwidths.quantity*np.sqrt(2*np.pi))
eintegrated = table.Column(name='Integrated Intensity error',
    data=((amplsK.quantity**2*evwidths.quantity**2) + (vwidths.quantity**2*eamplsK.quantity**2))**0.5)
jtok = table.Column(name='Jy/K', data=makecol('jtok', unit=u.Jy/u.K))
eu = table.Column(name='EU_K', data=makecol('EU_K', unit=u.Jy/u.K))
species = table.Column(name='Species', data=makecol('species'))
qn = table.Column(name='QNs', data=makecol('qns'))
deg = table.Column(name='deg', data=makecol('deg'))
Aij = table.Column(name='Aij', data=makecol('aul'))
flag = table.Column(name='Flag', data=makecol('flag'))
corner = table.Column(name='Corner', data=makecol('corner'))

vre = re.compile('v=([0-9]+)')
vstate = [int(vre.search(ss).groups()[0]) for ss in species]
jre = re.compile('J=([0-9]+)-([0-9]+)')
Ju,Jl = zip(*[map(int, jre.search(ss).groups()) for ss in species])
qnv = (Column(name='v', data=vstate))
qnju = (Column(name='J$_u$', data=Ju))
qnjl = (Column(name='J$_l$', data=Jl))


tbl1 = table.Table([linenames, species, qn, qnv, qnju, qnjl, freqs, velos,
    evelos, vwidths, evwidths, amplsK, eamplsK, integrated, eintegrated, jtok,
    eu, deg, Aij, flag, corner])

tbl1.write(paths.tpath('fitted_corner_lines.txt'), format='ascii.fixed_width', overwrite=True)



# create a subtable for the paper

linenames = table.Column([lines.texnames[ln]
                          if ln in lines.texnames
                          else ln
                          for ln, freq in zip(linenames, freqs)
                         ],
                         name='Line Name',
                        )


tbl = table.Table([linenames, qnv, qnju, qnjl, freqs, velos, evelos, vwidths, evwidths, amplsK, eamplsK, integrated, eintegrated, eu, flag, corner])

tbl.sort('Frequency')

bad_fits = []

badmask = np.array([ln in bad_fits for ln in linenames], dtype='bool')
badmask |= ((tbl['Fitted Width error'] > tbl['Fitted Width']) |
            (tbl['Fitted Velocity error'] > 5) |
            np.array([flg[1] in 'nq' for flg in tbl['Flag']])
           )

cornertbl = tbl
tbl.write(paths.tpath('corner_line_fits.txt'), format='ascii.fixed_width', overwrite=True)

formats = {'Frequency': lambda x: "{0:0.5f}".format(x),
           'Fitted Width': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Width error': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Velocity': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Velocity error': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Fitted Amplitude error K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Integrated Intensity': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'Integrated Intensity error': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'EU_K': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
          }
rename = {'Fitted Width':'Width',
          'Fitted Width error':'Width error',
          'Fitted Velocity':'Velocity',
          'Fitted Velocity error':'Velocity error',
          'Fitted Amplitude error K':'Amplitude error',
          'Fitted Amplitude K':'Amplitude',
          'Integrated Intensity': '$\int T_A dv$',
          'EU_K': 'E$_U$',
         }
maskcols = ['Fitted Width',
            'Fitted Width error',
            'Fitted Velocity',
            'Fitted Velocity error',
            'Fitted Amplitude error K',
            'Fitted Amplitude K',
            'Integrated Intensity',
            'Integrated Intensity error',
           ]

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
                if started_data:
                    continue
                else:
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
                fh.write("\hline\n")
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
                    vstate_lines.append(line)
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
#merge_errors(tbl, 'Integrated Intensity', 'Integrated Intensity error')
merge_errors(tbl, '$\int T_A dv$', 'Integrated Intensity error')
del formats['Velocity']
del formats['Width']
del formats['Amplitude']
del formats['$\int T_A dv$']

fulltbl = tbl
tbl = tbl['Line Name', 'v', 'J$_u$', 'J$_l$', 'Frequency', 'Velocity',
          'Amplitude', '$\int T_A dv$', 'E$_U$', 'Corner']


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
    latexdict['caption'] = 'Parameters of {0} lines obtained with Gaussian fits'.format(texsalt)
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
    #tbl[mask][columns].write(paths.texpath2('{0}_line_parameters.tex'.format(salt)),
    #                         formats=formats,
    #                         latexdict=latexdict,
    #                         overwrite=True)
    ##print(tbl[mask][columns])
    #assert len(tbl[mask][columns]) == mask.sum()
    #label_by_vstate(paths.texpath2('corner_{0}_line_parameters.tex'.format(salt)),
    #                ncols=len(columns)
    #               )

# abundance measurements

maskNaCl = np.array([ln.startswith(salt_to_barton['NaCl']) for ln in tbl['Line Name']])
maskNa37Cl = np.array([ln.startswith(salt_to_barton['Na37Cl']) for ln in tbl['Line Name']])
Na37Cltbl = tbl[maskNa37Cl]
NaCltbl = tbl[maskNaCl]

def get_meas(x):
    try:
        return float(x.split()[0])
    except:
        return np.nan

abund = {}
for row in NaCltbl:
    match = (Na37Cltbl['v'] == row['v']) & (Na37Cltbl['J$_u$'] == row['J$_u$']) & (Na37Cltbl['J$_l$'] == row['J$_l$'])
    if match.any():
        if match.sum() > 1:
            raise ValueError
        abund[row['Line Name']] = get_meas(row['Amplitude']) / get_meas(Na37Cltbl[match]['Amplitude'][0])


maskKCl = np.array([ln.startswith(salt_to_barton['KCl']) for ln in tbl['Line Name']])
maskK37Cl = np.array([ln.startswith(salt_to_barton['K37Cl']) for ln in tbl['Line Name']])
K37Cltbl = tbl[maskK37Cl]
KCltbl = tbl[maskKCl]

for row in KCltbl:
    match = (K37Cltbl['v'] == row['v']) & (K37Cltbl['J$_u$'] == row['J$_u$']) & (K37Cltbl['J$_l$'] == row['J$_l$'])
    if match.any():
        if match.sum() > 1:
            raise ValueError
        abund[row['Line Name']] = get_meas(row['Amplitude']) / get_meas(K37Cltbl[match]['Amplitude'][0])

print(abund)

print("Measured 35/37Cl abundance = {0} +/- {1}".format(np.mean(list(abund.values())),
                                                        np.std(list(abund.values()))))


maskK41Cl = np.array([ln.startswith(salt_to_barton['41KCl']) for ln in tbl['Line Name']])
K41Cltbl = tbl[maskK41Cl]
KCltbl = tbl[maskKCl]

for row in KCltbl:
    match = (K41Cltbl['v'] == row['v']) & (K41Cltbl['J$_u$'] == row['J$_u$']) & (K41Cltbl['J$_l$'] == row['J$_l$'])
    if match.any():
        if match.sum() > 1:
            raise ValueError
        abund[row['Line Name']+"_41"] = get_meas(row['Amplitude']) / get_meas(K41Cltbl[match]['Amplitude'][0])



# plot some things....

pl.figure(1).clf()

kclmask = np.array(['K-35Cl' in row['Species'] for row in tbl1])
k41clmask = np.array(['41K-35Cl' in row['Species'] for row in tbl1])
k37clmask = np.array(['K-37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][kclmask], tbl1['Fitted Amplitude K'][kclmask], 'o', label='KCl')
pl.plot(tbl1['EU_K'][k37clmask], tbl1['Fitted Amplitude K'][k37clmask], 's', label='K$^{37}$Cl')
pl.plot(tbl1['EU_K'][k41clmask], tbl1['Fitted Amplitude K'][k41clmask], 'd', label='$^{41}$KCl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [K]")
pl.legend(loc='best')
pl.savefig(paths.fpath('corners_KCl_amp_vs_eu.pdf'))

pl.figure(2).clf()
naclmask = np.array(['Na-35Cl' in row['Species'] for row in tbl1])
na37clmask = np.array(['Na-37Cl' in row['Species'] for row in tbl1])
pl.plot(tbl1['EU_K'][naclmask], tbl1['Fitted Amplitude K'][naclmask], 'o', label='NaCl')
pl.plot(tbl1['EU_K'][na37clmask], tbl1['Fitted Amplitude K'][na37clmask], 's', label='Na$^{37}$Cl')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Fitted amplitude [K]")
pl.legend(loc='best')
pl.savefig(paths.fpath('corners_NaCl_amp_vs_eu.pdf'))

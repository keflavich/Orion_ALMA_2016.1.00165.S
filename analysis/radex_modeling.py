import pyradex
import pylab as pl
import numpy as np
import paths
import dust_emissivity
from astropy import units as u
from astropy.table import Table, Column

rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
rslt = rr()

tbl = Table.read(paths.tpath('line_fits.txt'), format='ascii.fixed_width')
tbl.add_column(name='RADEX_row', col=Column(np.zeros(len(tbl), dtype='int')))

for row in tbl:
    upper = ("{0:<6}".format(f"{row['v']}_{row['J$_u$']}")).encode()
    lower = ("{0:<6}".format(f"{row['v']}_{row['J$_l$']}")).encode()
    match = (rslt['upperlevel'] == upper) & (rslt['lowerlevel'] == lower)
    if row['v'] <= 3:
        assert match.sum() == 1
    if match.sum() == 1:
        row['RADEX_row'] = np.argmax(match)
    else:
        row['RADEX_row'] = -999

fittbl = tbl[(tbl['RADEX_row'] >= 0) & (tbl['Fitted Amplitude error K'] <
                                        tbl['Fitted Amplitude K'])
            ]

def chi2(rdx, tbl=fittbl):
    inds = tbl['RADEX_row']

    # need non-background-subtracted brightness
    predicted_brightness = rdx['Tex'][inds] * (1-np.exp(-rdx['tau'][inds]))

    obs = tbl['Fitted Amplitude K']
    err = tbl['Fitted Amplitude error K']

    return (((obs-predicted_brightness)/err)**2).sum()

def GOF_plot(rdx, tbl=fittbl):

    inds = tbl['RADEX_row']
    predicted_brightness = rdx['Tex'][inds] * (1-np.exp(-rdx['tau'][inds]))

    obs = tbl['Fitted Amplitude K']
    err = tbl['Fitted Amplitude error K']

    energy = rdx['upperstateenergy'][inds]

    pl.plot(energy, predicted_brightness, 's')
    pl.errorbar(energy, obs, yerr=err, linestyle='none', marker='.')

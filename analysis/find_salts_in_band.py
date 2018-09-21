import json
import pyspeckit

from astropy import units as u
from astropy.utils.console import ProgressBar
from astropy.table import Table,Column

import paths
from salt_tables import salt_tables

salts_in_band = {}

for spw in (0,1,2,3):
    for band in ('B3', 'B6', 'B7.lb'):
        fn = paths.dpath('stacked_spectra/OrionSourceI_{band}_spw{0}_robust0.5.fits'
                         .format(spw, band=band))
        sp = pyspeckit.Spectrum(fn)

        for tbl in ProgressBar(salt_tables):
            for row in tbl:
                if sp.xarr.in_range(u.Quantity(row['Freq'], u.GHz)):
                    salts_in_band[row['Species']] = (row['Freq'], row['E_U'], band, spw)


with open('salts_in_band.json', 'w') as fh:
    json.dump(salts_in_band, fh)

tbl = Table([Column(name='Species', data=list(salts_in_band.keys())),
             Column(name='Frequency', data=[salts_in_band[key][0] for key in salts_in_band]),
             Column(name='E_U', data=[salts_in_band[key][1] for key in salts_in_band]),
             Column(name='Band', data=[salts_in_band[key][2] for key in salts_in_band]),
             Column(name='SPW', data=[salts_in_band[key][3] for key in salts_in_band]),
             Column(name='Flag', data=['----' for key in salts_in_band]),
            ])

tbl.sort('Frequency')
tbl.write(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')

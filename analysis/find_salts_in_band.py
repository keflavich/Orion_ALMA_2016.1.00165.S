import os
import json
import pyspeckit

from astropy import units as u
from astropy.utils.console import ProgressBar
from astropy.table import Table,Column

import paths
from salt_tables import salt_tables

salttablepath = paths.tpath('salts_in_band.ipac')

salts_in_band = {}

for spw in (0,1,2,3):
    for band in ('B3', 'B6', 'B7.lb'):
        fn = paths.dpath('stacked_spectra/OrionSourceI_{band}_spw{0}_robust0.5.fits'
                         .format(spw, band=band))
        sp = pyspeckit.Spectrum(fn)

        for tbl in ProgressBar(salt_tables):
            for row in tbl:
                if sp.xarr.in_range(u.Quantity(row['Freq'], u.GHz)):
                    salts_in_band[row['Species']] = (float(row['Freq']),
                                                     float(row['E_U']),
                                                     band,
                                                     spw,
                                                     float(row['Aij']))


with open('salts_in_band.json', 'w') as fh:
    json.dump(salts_in_band, fh)

tbl = Table([Column(name='Species', data=list(salts_in_band.keys())),
             Column(name='Frequency', data=[salts_in_band[key][0] for key in salts_in_band]),
             Column(name='E_U', data=[salts_in_band[key][1] for key in salts_in_band]),
             Column(name='Band', data=[salts_in_band[key][2] for key in salts_in_band]),
             Column(name='SPW', data=[salts_in_band[key][3] for key in salts_in_band]),
             Column(name='Aij', data=[salts_in_band[key][4] for key in salts_in_band]),
             Column(name='Flag', data=['----' for key in salts_in_band]),
            ])

tbl.sort('Frequency')
if not os.path.exists(salttablepath):
    tbl.write(salttablepath, format='ascii.ipac', overwrite=False)
else:
    tbl_ = Table.read(salttablepath, format='ascii.ipac')
    new_flagcolumn = [tbl_[tbl_['Species'] == row['Species']]['Flag'][0]
                      for row in tbl
                      if any(tbl_['Species'] == row['Species'])
                     ]
    tbl.remove_column('Flag')
    tbl.add_column(Column(name='Flag', data=new_flagcolumn))
    tbl.write(salttablepath, format='ascii.ipac', overwrite=True)

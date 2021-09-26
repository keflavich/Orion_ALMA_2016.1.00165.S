import paths
from astropy import units as u
from astropy.table import Table
import regions


tbl = Table.read(paths.rpath('3mmmSiOmasers_sourceI_average_centroids.txt'),
                 format='ascii', data_start=3)

center_reg = regions.Regions.read(paths.rpath('sourceI_center.reg'))[0]
center = center_reg.center[0]

with open(paths.rpath('3mm_sio_masers.reg'), 'w') as fh:
    fh.write('icrs\n')
    for row in tbl:
        ra = row[3]*u.arcsec + center.ra - 18.8*u.mas
        dec = row[5]*u.arcsec + center.dec + 116.8*u.mas
        fh.write('point({0}, {1}) # point=x\n'
                 .format(ra.to(u.deg).value,
                         dec.to(u.deg).value))

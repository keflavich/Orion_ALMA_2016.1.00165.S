"""
A dictionary of source names and their locations
"""

from astropy import coordinates, units as u

sources = {'SourceI': (83.81049333, -5.37517219),
           'Source8': (83.81072758, -5.375359983),
           'Source9': (83.81167742, -5.375188087),
           'Source1': (83.81106842, -5.377357087),
           'Source11': (83.80999342, -5.375112087),
           'Source18': (83.80939942, -5.374458087),
           'BN': (83.80878333, -5.37296361),
          }

def casafmt(crd):
    center = coordinates.SkyCoord(*crd, frame='icrs', unit=(u.deg, u.deg))

    return 'ICRS {} {}'.format(center.ra.to_string(unit=u.hour),
                               center.dec.to_string())

sources_fmtd = {key: casafmt(val) for key, val in sources.items()}

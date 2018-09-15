import os
import glob
import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import pyspeckit
import pylab as pl
from constants import vcen
from astroquery.splatalogue import Splatalogue
import lines

all_lines = {**lines.disk_lines, **lines.absorbers}

ided_linenames = sorted(all_lines.keys())
ided_linefreqs = u.Quantity([all_lines[x] for x in ided_linenames
                             if 'U' not in x])
ided_linetexnames = [lines.texnames[x] if x in lines.texnames else x
                     for x in ided_linenames
                     if 'U' not in x]

KCl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' KCl')
K37Cl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' K37Cl')
AlO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlO')
AlO = [row for row in AlO if len(row['Resolved QNs']) < 10]
NaCl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCl')
Na37Cl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' Na37Cl')

tables = [KCl, K37Cl, AlO, NaCl, Na37Cl]

def linename(row):
    return "{0} {1}".format(row['Species'], row['Resolved QNs'])
def freq(row):
    return u.Quantity(row['Freq-GHz(rest frame,redshifted)'] or
                      row['Meas Freq-GHz(rest frame,redshifted)'],
                      u.GHz)

linenames = [linename(row) for tbl in tables for row in tbl]
linetexnames = [linename(row) for tbl in tables for row in tbl]
linefreqs = [freq(row) for tbl in tables for row in tbl]

for fn in glob.glob(paths.dpath('stacked_spectra/OrionSourceI_*robust*fits')):
    print(fn)
    sp_st = pyspeckit.Spectrum(fn)
    sp_st.plotter(figure=pl.figure(0, figsize=(16,6)), clear=True)

    basefn = os.path.split(fn)[-1]

    sp_st.plotter(ymax=0.1)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/{0}'
                                      .format(basefn.replace("fits","pdf")))
                          )

    sp_st.plotter(ymin=-0.01, ymax=0.01)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)
    sp_st.plotter.line_ids(ided_linetexnames, ided_linefreqs, velocity_offset=-vcen,
                           plot_kwargs=dict(color='b'))
    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

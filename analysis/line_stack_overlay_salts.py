import numpy as np
import os
import glob
import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import pyspeckit
import pylab as pl
from constants import vcen
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
import lines

all_lines = {**lines.disk_lines, **lines.absorbers}

ided_linenames = sorted(all_lines.keys())
ided_linefreqs = u.Quantity([all_lines[x] for x in ided_linenames
                             #if 'U' not in x
                            ])
ided_linetexnames = [lines.texnames[x] if x in lines.texnames else x
                     for x in ided_linenames
                     #if 'U' not in x
                    ]

KCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' KCl'))
K37Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' K37Cl'))
K41Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='41KCl'))
AlCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlCl'))
#AlF = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlF'))
#NaF = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaF'))
#NaO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaO'))
#NaOH = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaOH'))
#NaCH = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCH'))
CaCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' CaCl'))
AlO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlO'))
#AlO = AlO[np.array([len(row['QNs']) < 10 for row in AlO])]
NaCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCl'))
Na37Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' Na37Cl'))
MgCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' MgCl'))
#MgCl = [row for row in MgCl if len(row['Resolved QNs']) < 20]
#not detected:
# HCl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' HCl')
#NaCN = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCN')
#NaO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaO')
#FeCO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' FeCO')
#whiffs at 230 GHz:
# MgCCH = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='MgCCH')
#too many isotopologues don't have anything
# SiS = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='Silicon monosulfide')
#LOTS of lines.  No way.
# MnO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='MnO')
# SO is real, but only the main isotopologue?
#SO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='Sulfur Monoxide')
SO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' SO '))
S34O = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' 34SO '))
# only low-J lines of SO2...
# no, SO2 isn't really believable.
SO2 = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' SO2',
                                 energy_max=500, energy_type='eu_k'))


salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl,]
tables = []

def linename(row):
    return "{0} {1}".format(row['Species'], row['QNs'])
def freq(row):
    return u.Quantity(row['Freq'], u.GHz)

linenames = [linename(row) for tbl in tables for row in tbl]
linetexnames = [linename(row) for tbl in tables for row in tbl] + ided_linetexnames
linefreqs = np.hstack([u.Quantity([freq(row) for tbl in tables for row in tbl], u.GHz).value,
                       ided_linefreqs.value])
linefreqs = u.Quantity(linefreqs, u.GHz)

flist = [fn] if 'fn' in locals() else glob.glob(paths.dpath('stacked_spectra/OrionSourceI_*robust0.5.fits'))
for fn in flist:
    print(fn)
    sp_st = pyspeckit.Spectrum(fn)
    pl.figure(0, figsize=(16,6)).clf()
    sp_st.plotter(figure=pl.figure(0, figsize=(16,6)), clear=True)

    basefn = os.path.split(fn)[-1]

    sp_st.plotter(ymax=0.1)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/{0}'
                                      .format(basefn.replace("fits","pdf")))
                          )

    sp_st.plotter(ymin=-0.005, ymax=0.01)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)
    #sp_st.plotter.line_ids(ided_linetexnames, ided_linefreqs, velocity_offset=-vcen,
    #                       plot_kwargs=dict(color='b'))
    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )
    for tbl in salt_tables:
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors='b', linestyles=':')

    for row in AlCl:
        frq = u.Quantity(row['Freq'], u.GHz).value
        if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
            sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                      -0.05, 0.10,
                                      colors='r', linestyles=':')

    for row in MgCl:
        frq = u.Quantity(row['Freq'], u.GHz).value
        if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
            sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                      -0.05, 0.10,
                                      colors='g', linestyles=':')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )


    for (a,b) in zip(linetexnames, linefreqs):
        if (b>sp_st.xarr.min()) and (b<sp_st.xarr.max()) and a not in ided_linetexnames:
            print("'{0}': {1}*u.{2},".format(a,b.value,b.unit))

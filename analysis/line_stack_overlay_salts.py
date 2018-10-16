import numpy as np
import os
import glob
import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import radio_beam
import pyspeckit
import pylab as pl
from constants import vcen
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
import lines
from salt_tables import (salt_tables, SO, SO2, HCl, sis_tables, AlCl, AlF, Al37Cl,
                         NaF, AlO, AlOH, NaCN)

all_lines = {**lines.disk_lines, **lines.absorbers}

ided_linenames = sorted(all_lines.keys())
ided_linefreqs = u.Quantity([all_lines[x] for x in ided_linenames
                             #if 'U' not in x
                            ])
ided_linetexnames = [lines.texnames[x] if x in lines.texnames else x
                     for x in ided_linenames
                     #if 'U' not in x
                    ]

#salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]
salt_colors = ['b', 'm', 'darkgreen', 'orange', 'c', 'y']
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

detection_table = table.Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')
nondetections = (detection_table['Flag'] == '-n') | (detection_table['Flag'] == 'cn')
detection_table = detection_table[~nondetections]

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

    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))
    #sp_st.plotter.line_ids(ided_linetexnames, ided_linefreqs, velocity_offset=-vcen,
    #                       plot_kwargs=dict(color='b'))
    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    #sp_st.plotter(ymin=-0.0025, ymax=0.01)
    # use the salt names directly.  This is for labeling of the colored
    # lines; the publication-ready stuff still uses lines.py
    #sp_st.plotter.line_ids(detection_table['Species'],
    #                       u.Quantity(detection_table['Frequency'], u.GHz),
    #                       velocity_offset=-vcen,
    #                       auto_yloc_fraction=0.8)

    for tbl,color in zip(salt_tables, salt_colors):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':')

    #for row in HCl:
    #    frq = u.Quantity(row['Freq'], u.GHz).value
    #    if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
    #        sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
    #                                  -0.05, 0.10,
    #                                  colors='g', linestyles='--')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    # Do another one just for SiO
    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))


    for tbl,color in zip(sis_tables, salt_colors):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/sis_diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    # Do another one just for alcl
    sp_st.plotter(ymin=-0.0025, ymax=0.01)

    # uses lines.py
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen,
                           auto_yloc_fraction=0.8)

    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))

    for tbl,color in zip([AlCl, AlF, Al37Cl, NaF, AlOH, AlO], salt_colors):
        for row in tbl:
            frq = u.Quantity(row['Freq'], u.GHz).value
            if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                          -0.05, 0.10,
                                          colors=color, linestyles=':')


    sp_st.plotter.savefig(paths.fpath('stacked_spectra/alcl_diagnostic_lines_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )


    for (a,b) in zip(linetexnames, linefreqs):
        if (b>sp_st.xarr.min()) and (b<sp_st.xarr.max()) and a not in ided_linetexnames:
            print("'{0}': {1}*u.{2},".format(a,b.value,b.unit))

    for speciesname, species in (('NaCN', NaCN), ('SO2',SO2), ('SO', SO)):
        # Do another one just for nacn
        sp_st.plotter(ymin=-0.0025, ymax=0.01)

        # uses lines.py
        sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen,
                               auto_yloc_fraction=0.8)

        for txt in sp_st.plotter.axis.texts:
            txt.set_backgroundcolor((1,1,1,0.9))

        for tbl,color in zip([species], ['b']):
            for row in tbl:
                frq = u.Quantity(row['Freq'], u.GHz).value
                if frq > sp_st.xarr.min().value and frq < sp_st.xarr.max().value:
                    sp_st.plotter.axis.vlines(frq*(1+vcen/constants.c).decompose().value,
                                              -0.05, 0.10,
                                              colors=color, linestyles=':')


        sp_st.plotter.savefig(paths.fpath('stacked_spectra/{1}_diagnostic_lines_labeled_{0}'
                                          .format(basefn.replace("fits","pdf"),
                                                  speciesname
                                                 ))
                             )

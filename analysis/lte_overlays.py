import numpy as np
import os
import glob
import paths
import radio_beam
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import pyspeckit
import pylab as pl
from constants import vcen
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
import lines
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_model
from parse_exomol_files import get_partition_func
from salt_tables import NaCl, Na37Cl, KCl, K37Cl, K41Cl, K41Cl37, salt_tables

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

flist = [fn] if 'fn' in locals() else glob.glob(paths.dpath('stacked_spectra/OrionSourceI_*robust0.5.fits'))
for fn in flist:

    basefn = os.path.split(fn)[-1]

    if 'B7' in basefn and 'lb' not in basefn:
        continue

    print(fn)

    beam = (radio_beam.Beam(0.1*u.arcsec, 0.08*u.arcsec) if 'B3' in fn else
            radio_beam.Beam(0.043*u.arcsec, 0.034*u.arcsec) if 'B6' in fn else
            radio_beam.Beam(0.029*u.arcsec, 0.022*u.arcsec))

    sp_st = pyspeckit.Spectrum(fn)

    jytok = beam.jtok(sp_st.xarr.mean())
    sp_st.data *= jytok.value
    sp_st.unit = u.K


    pl.figure(0, figsize=(16,6)).clf()
    sp_st.plotter(figure=pl.figure(0, figsize=(16,6)), clear=True,
                  ymin=-0.0025*jytok.value, ymax=0.01*jytok.value)



    for txt in sp_st.plotter.axis.texts:
        txt.set_backgroundcolor((1,1,1,0.9))


    for molname, molfullname, mol, col, tem, color in (
        ('NaCl', '23Na-35Cl', NaCl, 4e12, 1000, 'b'),
        ('NaCl', '23Na-37Cl', Na37Cl, 1e12, 1000, (0,0.1,0.9)),
        ('KCl', '39K-35Cl', KCl, 2e12, 1000, 'r'),
        ('KCl', '39K-37Cl', K37Cl, 0.5e12, 1000, (1,0.1,0)),
        ('KCl', '41K-35Cl', K41Cl, 0.5e12, 1000, (1,0.0,0.1)),
        ('KCl', '41K-37Cl', K41Cl37, 0.1e12, 1000, (1,0.1,0.1)),
       ):

        freqs = mol['Freq'].quantity

        match = (freqs > sp_st.xarr.min()) & (freqs < sp_st.xarr.max())
        freqs = freqs[match]
        aij = mol['Aij'][match].quantity.to(u.s**-1).value
        deg = mol['gu'][match]
        EU = mol['E_U'][match].quantity.to(u.erg, u.temperature_energy()).value

        def partfunc(tem):
            tbl = get_partition_func(molname, molfullname)
            return np.interp(tem, tbl['Temperature'], tbl['Q'])

        def mol_modfunc(xarr, vcen, width, tex, column):
            return generate_model(xarr, vcen, width, tex, column, freqs=freqs,
                                  aij=aij,
                                  deg=deg,
                                  EU=EU, partfunc=partfunc)

        model = mol_modfunc(sp_st.xarr,
                            -1*u.km/u.s,
                            4*u.km/u.s,
                            1000*u.K,
                            col*u.cm**-2)
        model[model<0.1] = np.nan

        sp_st.plotter.axis.plot(sp_st.xarr,
                                model,
                                linewidth=2,
                                alpha=0.75,
                                zorder=-5,
                                color=color,
                               )

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lte_overlay_K_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen,
                           auto_yloc_fraction=0.8)

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lte_overlay_K_labeled_{0}'
                                      .format(basefn.replace("fits","pdf")))
                         )

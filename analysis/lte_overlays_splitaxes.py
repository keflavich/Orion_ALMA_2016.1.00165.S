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
for band in (3,6,7):
    fig = pl.figure(0, figsize=(16,6))


    flist_b = [fn for fn in flist
               if 'B{0}'.format(band) in fn
               and ('.lb' in fn if 'B7' in fn else True)
              ]

    for flist_c,suffix in ((flist_b[:2], 'a'), (flist_b[2:], 'b')):
        fig.clf()
        fig, axes = pl.subplots(1, 2, sharey=True, num=0)

        for fn,ax in zip(flist_c, axes):

            ax.cla()

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

            sp_st.baseline.includemask = np.abs(sp_st.data) < 10
            sp_st.baseline(order=3)

            ax.plot(sp_st.xarr, sp_st.data, linewidth=0.5, color='k',
                    linestyle='steps-mid')
            #ax.set_xlabel("Frequency (GHz)")
            ax.set_ylabel("Brightness Temperature (K)")


            for txt in ax.texts:
                txt.set_backgroundcolor((1,1,1,0.9))


            for molname, molfullname, mol, col, tem, color in (
                ('NaCl', '23Na-35Cl', NaCl, 4e12, 1000, 'g'),
                ('NaCl', '23Na-37Cl', Na37Cl, 2e12, 1000, (0,0.9,0.1)),
                ('KCl', '39K-35Cl', KCl, 2e12, 1000, 'r'),
                ('KCl', '39K-37Cl', K37Cl, 1e12, 1000, (1,0.1,0)),
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
                                    -vcen,
                                    4*u.km/u.s,
                                    1000*u.K,
                                    col*u.cm**-2)
                model[model==0] = np.nan

                ax.plot(sp_st.xarr,
                        model,
                        linewidth=1,
                        alpha=0.75,
                        zorder=5,
                        color=color,
                       )
            ax.set_ylim(-20, 60)

        axes[0].spines['right'].set_visible(False)
        #axes[1].spines['left'].set_visible(False)
        #axes[1].spines['right'].set_visible(False)
        #axes[2].spines['left'].set_visible(False)
        #axes[2].spines['right'].set_visible(False)
        axes[1].spines['left'].set_visible(False)
        axes[0].yaxis.tick_left()
        #axes[1].yaxis.set_ticks([])
        #axes[2].yaxis.set_ticks([])
        axes[1].yaxis.tick_right()
        axes[1].set_ylabel('')
        #axes[2].set_ylabel('')
        #axes[3].set_ylabel('')
        pl.subplots_adjust(wspace=0.04)
        pl.suptitle("Frequency (GHz)", y=0.05)

        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict( color='k', clip_on=False)
        axes[0].plot((1 - d, 1 + d), (1-d, 1+d), transform=axes[0].transAxes, **kwargs)  # top-right diagonal
        axes[0].plot((1 - d, 1 + d), (-d, +d), transform=axes[0].transAxes, **kwargs)  # bottom-right diagonal
        #axes[1].plot((1 - d, 1 + d), (1-d, 1+d), transform=axes[1].transAxes, **kwargs)  # top-right diagonal
        #axes[1].plot((1 - d, 1 + d), (-d, +d), transform=axes[1].transAxes, **kwargs)  # bottom-right diagonal
        axes[1].plot(( - d,  + d), (1-d, 1+d), transform=axes[1].transAxes, **kwargs)  # top-left diagonal
        axes[1].plot(( - d,  + d), (-d, +d), transform=axes[1].transAxes, **kwargs)  # bottom-left diagonal
        #axes[2].plot((1 - d, 1 + d), (1-d, 1+d), transform=axes[2].transAxes, **kwargs)  # top-right diagonal
        #axes[2].plot((1 - d, 1 + d), (-d, +d), transform=axes[2].transAxes, **kwargs)  # bottom-right diagonal
        #axes[2].plot(( - d,  + d), (1-d, 1+d), transform=axes[2].transAxes, **kwargs)  # top-left diagonal
        #axes[2].plot(( - d,  + d), (-d, +d), transform=axes[2].transAxes, **kwargs)  # bottom-left diagonal
        #axes[3].plot(( - d,  + d), (1-d, 1+d), transform=axes[3].transAxes, **kwargs)  # top-left diagonal
        #axes[3].plot(( - d,  + d), (-d, +d), transform=axes[3].transAxes, **kwargs)  # bottom-left diagonal

        #kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        #ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        #ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

        pl.savefig(paths.fpath('stacked_spectra/squash_lte_overlay_K_{1}_{0}'
                               .format(basefn.replace("fits","pdf"), suffix)),
                   bbox_inches='tight'
                  )

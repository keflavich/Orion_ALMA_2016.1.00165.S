import paths
import lines
import imp
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import pyspeckit
import pylab as pl
from constants import vcen
imp.reload(lines)

all_lines = {**lines.disk_lines, **lines.absorbers}

linenames = sorted(all_lines.keys())
linefreqs = u.Quantity([all_lines[x] for x in linenames])
linetexnames = [lines.texnames[x] for x in linenames]

for spw in (0,1,2,3):
    sp_st = pyspeckit.Spectrum(paths.dpath('stacked_spectra/OrionSourceI_B6_spw{0}.fits'.format(spw)))
    sp_hc = pyspeckit.Spectrum(paths.dpath('stacked_spectra/hotcore_spectrum_full_OrionSourceI_B6_maskedclean_spw{0}_lines.fits'.format(spw)))
    sp_hc.xarr.convert_to_unit(u.GHz)
    sp_st.plotter(figure=pl.figure(spw, figsize=(16,6)), clear=True)
    sp_hc.plotter(axis=sp_st.plotter.axis, clear=False, color='b')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/hotcore_overlay_spw{0}_B6.pdf'.format(spw)))

    sp_st.plotter()
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_spw{0}_B6.pdf'.format(spw)))

    sp_st.plotter(ymin=-0.01, ymax=0.01)
    sp_st.plotter.line_ids(linetexnames, linefreqs, velocity_offset=-vcen)
    sp_st.plotter.savefig(paths.fpath('stacked_spectra/lines_labeled_spw{0}_B6_yzoom.pdf'.format(spw)))

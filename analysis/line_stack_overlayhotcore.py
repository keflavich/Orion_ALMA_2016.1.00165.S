import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import pyspeckit


for spw in (0,1,2,3):
    sp_st = pyspeckit.Spectrum(paths.dpath('stacked_spectra/OrionSourceI_B6_spw{0}.fits'.format(spw)))
    sp_hc = pyspeckit.Spectrum(paths.dpath('stacked_spectra/hotcore_spectrum_full_OrionSourceI_B6_maskedclean_spw{0}_lines.fits'.format(spw)))
    sp_hc.xarr.convert_to_unit(u.GHz)
    sp_st.plotter()
    sp_hc.plotter(axis=sp_st.plotter.axis, clear=False, color='b')

    sp_st.plotter.savefig(paths.fpath('stacked_spectra/hotcore_overlay_spw{0}_B6.pdf'.format(spw)))

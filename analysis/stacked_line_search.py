import numpy as np
import os
import spectral_cube.analysis_utilities
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.io import fits
import paths
from paths import fcp
import pylab as pl
import regions
import reproject

# step 1: create a velocity map

vmap_name = 'disk_velocity_map.fits'
if not os.path.exists(vmap_name):
    cube = SpectralCube.read(paths.dpath('cubes/OrionSourceI_Unknown_4_robust0.5.maskedclarkclean10000_medsub_K.fits'))
    cube = SpectralCube.read(paths.dpath('cubes/OrionSourceI_Unknown_4_robust0.5maskedclarkclean10000_medsub_K.fits'))
    m1 = cube.moment1()
    m0 = cube.moment0()
    mask = m0.value > 300

    vmap = m1
    vmap[~mask] = np.nan

    r =regions.read_ds9(paths.rpath('sourceI_enclosing_ellipse.reg'))[0]
    rp = r.to_pixel(vmap.wcs)
    mask = rp.to_mask()

    vmap_ = np.empty(vmap.shape)*np.nan
    vmap_[mask.bbox.slices] = vmap[mask.bbox.slices].value * mask.data
    hdu = vmap.hdu
    hdu.data = vmap_
    hdu.writeto(vmap_name, overwrite=True)
else:
    hdu = fits.open(vmap_name)[0]
vmap = spectral_cube.lower_dimensional_structures.Projection.from_hdu(hdu)


# step 2: stack

for band in ('B3', 'B6', 'B7'):
    for spw in (0,1,2,3):
        for robust in (-2, 0.5, 2):

            suffix = '.lb' if band == 'B7' else ''

            try:
                fn = fcp('OrionSourceI_only.{1}{3}.robust{2}.spw{0}.maskedclarkclean10000_medsub.image.pbcor.fits'
                         .format(spw, band, robust, suffix))
                fullcube = (SpectralCube.read(fn))
            except FileNotFoundError:
                fn = fcp('OrionSourceI_only.{1}{3}.robust{2}.spw{0}.clarkclean10000_medsub.image.pbcor.fits'
                         .format(spw, band, robust, suffix))
                fullcube = (SpectralCube.read(fn))
            print(fn,fullcube.spectral_extrema)
            fullcube = fullcube.with_spectral_unit(u.km/u.s,
                                                   velocity_convention='radio',
                                                   rest_value=fullcube.spectral_axis.mean())

            vmap_proj,_ = reproject.reproject_interp(vmap.hdu,
                                                     fullcube.wcs.celestial,
                                                     shape_out=fullcube.shape[1:])
            vmap_proj = u.Quantity(vmap_proj, u.km/u.s)

            stack = spectral_cube.analysis_utilities.stack_spectra(fullcube, vmap_proj,
                                                                   v0=0.0*u.km/u.s)
            fstack = stack.with_spectral_unit(u.GHz)

            fstack.write(paths.dpath('stacked_spectra/OrionSourceI_{1}{3}_spw{0}_robust{2}.fits'
                                     .format(spw, band, robust, suffix)),
                         overwrite=True)

            pl.clf()
            fstack.quicklook(filename=paths.fpath('stacked_spectra/OrionSourceI_{1}{3}_spw{0}_robust{2}.pdf')
                             .format(spw, band, robust, suffix))

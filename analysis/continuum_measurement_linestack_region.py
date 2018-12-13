import numpy as np
from spectral_cube import SpectralCube, wcs_utils, tests, Projection
import radio_beam
from astropy import convolution
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy import stats
import paths
import pylab as pl
from scipy.ndimage import map_coordinates
import scipy.signal
import reproject

from files import b3_hires_cont, b6_hires_cont, b7_hires_cont
from constants import source, extraction_path, origin, central_freqs


# vmap produced by stacked_line_search.py
vmap_name = paths.dpath('disk_velocity_map.fits')
hdu = fits.open(vmap_name)[0]
vmap = Projection.from_hdu(hdu)

b3beam = radio_beam.Beam.from_fits_header(fits.getheader(paths.dpath(b3_hires_cont)))

print("per-band continuum measurements in the spectral extraction aperture: ")
for ii,contfn in enumerate((b3_hires_cont, b6_hires_cont, b7_hires_cont)):
    band = contfn[14:16]

    conthdu = fits.open(paths.dpath(contfn))[0]

    ww = wcs.WCS(conthdu.header)

    #vmap_proj,_ = reproject.reproject_interp(vmap.hdu,
    #                                         ww,
    #                                         shape_out=conthdu.data.shape)
    #vmap_proj = u.Quantity(vmap_proj, u.km/u.s)
    #vmap_mask = np.isfinite(vmap_proj)
    cont_proj,_ = reproject.reproject_interp(conthdu,
                                             vmap.wcs,
                                             shape_out=vmap.shape)
    vmap_mask = np.isfinite(vmap)

    beam = radio_beam.Beam.from_fits_header(conthdu.header)
    jtok = beam.jtok(central_freqs[band])

    flux = cont_proj[vmap_mask].sum()
    intensity = np.nanmean(cont_proj[vmap_mask]) * jtok

    print(f"{band}: f={flux} Jy, <T_B>={intensity}")

    try:
        conv_beam = b3beam.deconvolve(beam)
    except ValueError:
        # b3beam can't deconvolve itself
        continue

    pixscale = wcs.utils.proj_plane_pixel_scales(vmap.wcs)[0]*u.deg
    kernel = conv_beam.as_kernel(pixscale)
    cont_proj_smooth = convolution.convolve_fft(cont_proj*jtok, kernel, preserve_nan=True)
    jtokb3 = b3beam.jtok(central_freqs[band])

    flux = (cont_proj_smooth/jtokb3)[vmap_mask].sum()
    intensity = np.nanmean(cont_proj_smooth[vmap_mask])

    print(f"Convolved to B3: {band}: f={flux} Jy, <T_B>={intensity}")

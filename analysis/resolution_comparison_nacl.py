import numpy as np
import pylab as pl
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import convolution
import radio_beam
import reproject

from files import b6_hires_cont, b7_hires_cont, b3_hires_cont
from constants import central_freqs
import paths

b3beam = radio_beam.Beam.from_fits_header(fits.getheader(paths.dpath(b3_hires_cont)))

fig = pl.figure(1, figsize=(10,10))
pl.clf()

#\includegraphics[scale=1,width=2.25in]{figures/OrionSourceI_NaClv=3_7-6_robust0.5.maskedclarkclean10000_medsub_K_peak_offset_contours.pdf}
#\includegraphics[scale=1,width=2.25in]{figures/OrionSourceI_NaClv=1_18-17_robust0.5.maskedclarkclean10000_medsub_K_peak_offset_contours.pdf}
#\includegraphics[scale=1,width=2.25in]{figures/OrionSourceI_NaClv=2_26-25_robust0.5.maskedclarkclean10000_medsub_K_peak_offset_contours.pdf}

sourcename = 'SourceI'
for ii,(linename, robust, band) in enumerate((('NaClv=3_7-6', 0.5, 'B3',),
                                             ('NaClv=1_18-17', 0.5, 'B6',),
                                             ('NaClv=2_26-25', 0.5, 'B7'),)):

    mx_fn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.fits').format(linename, sourcename, robust=robust)

    fh = fits.open(mx_fn)

    ww = wcs.WCS(fh[0].header)
    beam = radio_beam.Beam.from_fits_header(fh[0].header)

    jtok = beam.jtok(central_freqs[band])

    ax = pl.subplot(3, 3, ii+1)
    im = ax.imshow(fh[0].data-np.percentile(fh[0].data, 10), interpolation='none', origin='lower', cmap='viridis', vmin=0, vmax=300)
    ax.set_title(linename.replace("_"," "))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    if band == 'B7':
        cbax = fig.add_axes([0.99,0.655,0.02,0.3])
        cb = pl.colorbar(mappable=im, cax=cbax)
        cb.set_label("T$_B$ [K]")


    try:
        conv_beam = b3beam.deconvolve(beam)
    except ValueError:
        b3fn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.fits').format(linename, sourcename, robust=robust)
        b3fh = fits.open(b3fn)
        ax2 = pl.subplot(3, 3, ii+4)
        ax2.imshow(b3fh[0].data - np.percentile(b3fh[0].data, 10), interpolation='none', origin='lower', cmap='viridis', vmin=0, vmax=175)
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        # b3beam can't deconvolve itself
        continue

    ax2 = pl.subplot(3, 3, ii+4)

    pixscale = wcs.utils.proj_plane_pixel_scales(ww)[0]*u.deg
    kernel = conv_beam.as_kernel(pixscale)
    smoothed = convolution.convolve_fft(fh[0].data, kernel, preserve_nan=True)

    im = ax2.imshow(smoothed-np.percentile(smoothed, 10), interpolation='none', origin='lower', cmap='viridis', vmin=0, vmax=175)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    if band == 'B7':
        cbax = fig.add_axes([0.99,0.34,0.02,0.3])
        cb = pl.colorbar(mappable=im, cax=cbax)
        cb.set_label("T$_B$ [K]")

    ax3 = pl.subplot(3, 3, ii+7)

    reproj,_ = reproject.reproject_interp((smoothed-np.percentile(smoothed, 10), fh[0].header), b3fh[0].header)

    im = ax3.imshow(b3fh[0].data - np.percentile(b3fh[0].data, 10) - reproj, interpolation='none', origin='lower', cmap='viridis', vmin=-75, vmax=75)
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])

    if band == 'B7':
        cbax = fig.add_axes([0.99,0.025,0.02,0.3])
        cb = pl.colorbar(mappable=im, cax=cbax)
        cb.set_label("T$_B$ [K]")

pl.tight_layout()
pl.subplots_adjust(hspace=0.02, wspace=0.02)

pl.savefig(paths.fpath("resolution_comparison_nacl_peaks.pdf"), bbox_inches='tight')

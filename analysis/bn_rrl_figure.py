import os
import numpy as np
import spectral_cube
from spectral_cube import SpectralCube
import radio_beam
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.convolution import Gaussian1DKernel
from astropy.io import fits
import paths
import pylab as pl

# just to be sure...
pl.ion()
pl.close('all')

cube_B6 = SpectralCube.read(paths.dpath('cubes/OrionBN_only.B6.robust0.5.spw1.maskedclarkclean10000.image.pbcor.fits'))
rrl30a = cube_B6.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=231.900928*u.GHz).spectral_slab(-50*u.km/u.s, 100*u.km/u.s)
cube_B3 = SpectralCube.read(paths.dpath('cubes/OrionBN_only.B3.robust0.5.spw2.clarkclean10000.image.pbcor.fits'))
rrl40a = cube_B3.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=99.022952*u.GHz).spectral_slab(-50*u.km/u.s, 100*u.km/u.s)


h30a_cont_fn = paths.dpath('cubes/OrionBN_B6_robust0.5_spw1_cont.fits')
if os.path.exists(h30a_cont_fn):
    cube_cont_B6 = spectral_cube.lower_dimensional_structures.Projection.from_hdu(fits.open(h30a_cont_fn))
else:
    cube_cont_B6 = cube_B6.mask_out_bad_beams(0.1).median(axis=0)
    cube_cont_B6.write(h30a_cont_fn)

h40a_cont_fn = paths.dpath('cubes/OrionBN_B3_robust0.5_spw2_cont.fits')
if os.path.exists(h40a_cont_fn):
    cube_cont_B3 = spectral_cube.lower_dimensional_structures.Projection.from_hdu(fits.open(h40a_cont_fn))
else:
    cube_cont_B3 = cube_B3.mask_out_bad_beams(0.1).median(axis=0)
    cube_cont_B3.write(h40a_cont_fn)

h30a_contsub_fn = paths.dpath('cubes/OrionBN_H30a_contsub.fits')
if os.path.exists(h30a_contsub_fn):
    contsub_rrl30a = SpectralCube.read(h30a_contsub_fn)
else:


    contsub_rrl30a = rrl30a-cube_cont_B6

    contsub_rrl30a.write(h30a_contsub_fn)

h40a_contsub_fn = paths.dpath('cubes/OrionBN_H40a_contsub.fits')
if os.path.exists(h40a_contsub_fn):
    contsub_rrl40a = SpectralCube.read(h40a_contsub_fn)
else:


    contsub_rrl40a = rrl40a-cube_cont_B3

    contsub_rrl40a.write(h40a_contsub_fn)


ds_rrl30a_fn = paths.dpath('cubes/OrionBN_H30a_downsampled5kms_contsub.fits')
if not os.path.exists(ds_rrl30a_fn):
    # this is my by-hand estimated reasonable common beam
    onebeamcube_rrl30a = contsub_rrl30a.mask_out_bad_beams(0.1).convolve_to(radio_beam.Beam(0.059*u.arcsec, 0.047*u.arcsec, -86.0*u.deg))
    # 1.26 ~ channel width
    sm_contsub_rrl30a = onebeamcube_rrl30a.spectral_smooth(Gaussian1DKernel(((5/2.35)**2 - 1.26**2)**0.5))

    downsampled_contsub_rrl30a = sm_contsub_rrl30a.spectral_interpolate(np.arange(-25, 85, 5)*u.km/u.s)

    downsampled_contsub_rrl30a.write(ds_rrl30a_fn)
else:
    downsampled_contsub_rrl30a = SpectralCube.read(ds_rrl30a_fn)


ds_rrl40a_fn = paths.dpath('cubes/OrionBN_H40a_downsampled5kms_contsub.fits')
if not os.path.exists(ds_rrl40a_fn):
    # this is my by-hand estimated reasonable common beam
    onebeamcube_rrl40a = contsub_rrl40a.mask_out_bad_beams(0.1).convolve_to(radio_beam.Beam(0.090*u.arcsec, 0.061*u.arcsec, 43.6*u.deg))

    # already smoothed enough for this
    downsampled_contsub_rrl40a = onebeamcube_rrl40a.spectral_interpolate(np.arange(-25, 85, 5)*u.km/u.s)

    downsampled_contsub_rrl40a.write(ds_rrl40a_fn)
else:
    downsampled_contsub_rrl40a = SpectralCube.read(ds_rrl40a_fn)



mask = ((downsampled_contsub_rrl30a.spectral_axis > -10*u.km/u.s) & (downsampled_contsub_rrl30a.spectral_axis < 50*u.km/u.s))
whmask = np.where(mask)[0]

cont_B6 = spectral_cube.lower_dimensional_structures.Projection.from_hdu(fits.open(paths.dpath('Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5.image.tt0.pbcor.fits')))

BNcutout_cont_B6 = cont_B6[5513:5618,5037:5154]
BNcutout_cont_B6.quicklook()
BNcutout_cont_B6.FITSFigure.show_grayscale(stretch='linear')

cube_cont_B6.quicklook()
FFB6 = cube_cont_B6.FITSFigure

#for tx in FFB6.ax.texts:
#    tx.set_visible(False)
#
#for lyr in list(FFB6._layers):
#    FFB6.remove_layer(lyr)

for ii,ind in enumerate(whmask):
    ind = int(ind)
    vel = downsampled_contsub_rrl30a.spectral_axis[ind]
    slc = downsampled_contsub_rrl30a[ind,:,:]
    color = pl.cm.Spectral((ind-whmask.min())/mask.sum())
    FFB6.show_contour(slc.hdu, levels=[slc.max().value*0.95], colors=[color])
    FFB6.ax.annotate("{0:0.1f}".format(vel),
                     (0.20, 0.88-0.02*ii),
                     xycoords='figure fraction',
                     color=color)

#FFB6.ax.axis((41.559800470110005, 76.83040770863137, 38.045746948239589, 63.824569975782737))
FFB6.recenter(83.80877507, -5.372957363, radius=(0.10*u.arcsec).to(u.deg).value)
FFB6.save(paths.fpath('rrls/BN_H30a_velocity_contours_on_BNcontinuum.pdf'), dpi=150)
FFB6.save(paths.fpath('rrls/BN_H30a_velocity_contours_on_BNcontinuum.png'))







mask = ((downsampled_contsub_rrl40a.spectral_axis > -10*u.km/u.s) & (downsampled_contsub_rrl40a.spectral_axis < 50*u.km/u.s))
whmask = np.where(mask)[0]

cont_B3 = spectral_cube.lower_dimensional_structures.Projection.from_hdu(fits.open(paths.dpath('Orion_SourceI_B3_continuum_r-2.clean0.1mJy.image.tt0.pbcor.fits')))

BNcutout_cont_B3 = cont_B3[5569:5627,5333:5397]
BNcutout_cont_B3.quicklook()
BNcutout_cont_B3.FITSFigure.show_grayscale(stretch='linear')

cube_cont_B3.quicklook()
FFB3 = cube_cont_B3.FITSFigure

#for tx in FFB3.ax.texts:
#    tx.set_visible(False)
#
#for lyr in list(FFB3._layers):
#    FFB3.remove_layer(lyr)

for ii,ind in enumerate(whmask):
    ind = int(ind)
    vel = downsampled_contsub_rrl40a.spectral_axis[ind]
    slc = downsampled_contsub_rrl40a[ind,:,:]
    color = pl.cm.Spectral((ind-whmask.min())/mask.sum())
    FFB3.show_contour(slc.hdu, levels=[slc.max().value*0.95], colors=[color])
    FFB3.ax.annotate("{0:0.1f}".format(vel),
                     (0.20, 0.88-0.02*ii),
                     xycoords='figure fraction',
                     color=color)

#FFB3.ax.axis((41.559800470110005, 76.84040770863137, 38.045746948239589, 63.824569975782737))
FFB3.recenter(83.80877507, -5.372957363, radius=(0.10*u.arcsec).to(u.deg).value)
FFB3.save(paths.fpath('rrls/BN_H40a_velocity_contours_on_BNcontinuum.pdf'), dpi=150)
FFB3.save(paths.fpath('rrls/BN_H40a_velocity_contours_on_BNcontinuum.png'))




for tx in BNcutout_cont_B3.FITSFigure.ax.texts:
    tx.set_visible(False)

for lyr in list(BNcutout_cont_B3.FITSFigure._layers):
    BNcutout_cont_B3.FITSFigure.remove_layer(lyr)

BNcutout_cont_B3.FITSFigure.show_contour(cube_cont_B3.hdu, levels=np.array([0.2,0.5,0.8,0.9,0.95]) * cube_cont_B3.max().value, colors=['r']*6)


for tx in BNcutout_cont_B6.FITSFigure.ax.texts:
    tx.set_visible(False)

for lyr in list(BNcutout_cont_B6.FITSFigure._layers):
    BNcutout_cont_B6.FITSFigure.remove_layer(lyr)

BNcutout_cont_B6.FITSFigure.show_contour(cube_cont_B6.hdu, levels=np.array([0.2,0.5,0.8,0.9,0.95]) * cube_cont_B6.max().value, colors=['r']*6)

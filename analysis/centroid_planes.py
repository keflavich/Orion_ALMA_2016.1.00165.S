import numpy as np
import regions
from astropy import units as u
from astropy.modeling import models, fitting
from astropy import wcs
from spectral_cube import SpectralCube
from gaussfit_catalog.core import gaussfit_image
import paths
import os
import pylab as pl


regs = regions.read_ds9('../regions/velo_centroid_guesses_Unknown_4.reg')

guesses = {}

for reg in regs:
    vel = float(reg.meta['text'].strip("{}"))
    if vel in guesses:
        guesses[vel].append(reg)
    else:
        guesses[vel] = [reg]

velocities = np.array(sorted(guesses.keys()))*u.km/u.s

cubefn = '../FITS/cubes/OrionSourceI_Unknown_4_robust0.5.maskedclarkclean10000_medsub_K.fits'
basename = os.path.splitext(os.path.basename(cubefn))[0]
cube = SpectralCube.read(cubefn)
vdiff = np.abs(np.diff(cube.spectral_axis).mean())
rms = cube.std()
weights = np.ones(cube.shape[1:]) * rms.value**2

try:
    beam = cube.beam
except AttributeError:
    beam = cube.average_beams(1)
pixscale = wcs.utils.proj_plane_pixel_area(cube.wcs.celestial)**0.5 * u.deg
bmmaj_px = (beam.major / pixscale).decompose()
bmmin_px = (beam.minor / pixscale).decompose()
max_radius_in_beams = 1.25
max_offset_in_beams = 1.5

STDDEV_TO_FWHM = np.sqrt(8*np.log(2))

results = {}

for vel,vslice in zip(cube.spectral_axis, cube):
    closest = np.argmin(np.abs(vel-velocities))
    if np.abs(velocities[closest] - vel) > vdiff:
        #print("Skipping velocity {0}, closest is {1} -> {2}".format(vel, closest,
        #                                                            velocities[closest]))
        continue

    thisvel = velocities[closest].value
    guess_regs = guesses[thisvel]

    ampguess = vslice.max().value

    model_list = []
    for reg in guess_regs:

        p_init = models.Gaussian2D(amplitude=ampguess,
                                   x_mean=reg.center.x,
                                   y_mean=reg.center.y,
                                   x_stddev=bmmaj_px/STDDEV_TO_FWHM*0.75,
                                   y_stddev=bmmin_px/STDDEV_TO_FWHM*0.75,
                                   theta=beam.pa,
                                   bounds={'x_stddev':(bmmin_px/STDDEV_TO_FWHM*0.5,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'y_stddev':(bmmin_px/STDDEV_TO_FWHM*0.5,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'x_mean':(reg.center.x-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                     reg.center.x+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                           'y_mean':(reg.center.y-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                     reg.center.y+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                           'amplitude':(ampguess*0.1, ampguess*2.1)
                                          }
                                  )
        model_list.append(p_init)

    composite_model = model_list[0]
    if len(model_list) > 1:
        for mod in model_list[1:]:
            composite_model += mod

    fit_result = gaussfit_image(vslice.value, composite_model,
                                weights=weights, plot=True)
    results[thisvel] = fit_result

    if not os.path.exists(paths.fpath('velcentroid/diagnostics/{0}/'.format(basename))):
        os.mkdir(paths.fpath('velcentroid/diagnostics/{0}/').format(basename))
    pl.savefig(paths.fpath('velcentroid/diagnostics/{0}/{1}.png')
               .format(basename, thisvel))

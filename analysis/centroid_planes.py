import numpy as np
import regions
from astropy import units as u
from astropy.modeling import models, fitting
from astropy import wcs
from astropy import coordinates
from spectral_cube import SpectralCube
from gaussfit_catalog.core import gaussfit_image
import paths
import os
import pylab as pl
import shapely.geometry as geom


regs = regions.read_ds9(paths.rpath('velo_centroid_guesses_Unknown_4.reg'))

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


diskend_regs = regions.read_ds9(paths.rpath('diskends.reg'))
diskends = coordinates.SkyCoord([reg.center for reg in diskend_regs])

center = coordinates.SkyCoord(diskends.ra.mean(), diskends.dec.mean(), frame=diskends.frame)
ref_cen_x, ref_cen_y = cube.wcs.celestial.wcs_world2pix(center.ra.deg, center.dec.deg, 0)

def offset_to_point(xx, yy):
    line = geom.LineString(zip(diskends.ra.deg, diskends.dec.deg))
    point = geom.Point(xx, yy)
    return line.project(point)

ref_offset = offset_to_point(center.ra.deg, center.dec.deg)


fig1 = pl.figure(1)
pl.clf()
ax1 = fig1.add_subplot(1,2,1)
ax2 = fig1.add_subplot(1,2,2)

diskends_x, diskends_y = cube.wcs.celestial.wcs_world2pix(diskends.ra.deg,
                                                          diskends.dec.deg,
                                                          0)
ax1.plot(diskends_x-ref_cen_x, diskends_y-ref_cen_y, 'k-', linewidth=3, alpha=0.2, zorder=-100)
ax1.plot(0, 0, 'ko', markersize=5, alpha=0.5, zorder=-50)

cmap = pl.cm.Spectral_r
cmap = pl.cm.spectral
vmin, vmax = -15, 27
norm = pl.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

for vel in results:
    fit_result = results[vel]
    fitted = fit_result[0]
    n_submodels = fitted.n_submodels()
    color = cmap(norm(vel))
    if n_submodels > 1:
        for ii in range(n_submodels):
            xcen, ycen = getattr(fitted, 'x_mean_{0}'.format(ii)), getattr(fitted, 'y_mean_{0}'.format(ii))
            ax1.plot(xcen-ref_cen_x,
                     ycen-ref_cen_y,
                     color=color,
                     marker='o')
            ra,dec = cube.wcs.celestial.wcs_pix2world(xcen, ycen, 0)
            offset = offset_to_point(ra, dec)
            ax2.plot((offset-ref_offset)*3600, vel, color=color, marker='s')
    else:
        xcen, ycen = fitted.x_mean, fitted.y_mean
        ax1.plot(xcen-ref_cen_x, ycen-ref_cen_y, color=color, marker='o')
        ra,dec = cube.wcs.celestial.wcs_pix2world(xcen, ycen, 0)
        offset = offset_to_point(ra, dec)
        ax2.plot((offset-ref_offset)*3600, vel, color=color, marker='s')

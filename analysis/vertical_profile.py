import numpy as np
from spectral_cube import SpectralCube, wcs_utils, tests, Projection
import radio_beam
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

sourcename='SourceI'
robust = 0.5
linename='NaClv=2_26-25'
linename='NaClv=1_18-17' # use this to compute the continuum in the spectral extraction region
naclfn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.fits').format(linename, sourcename, robust=robust)
nacldata = fits.getdata(naclfn)
ww = wcs.WCS(fits.getheader(naclfn))
pixscale = wcs.utils.proj_plane_pixel_scales(ww)[0]

yy,xx = np.indices(np.array(nacldata.shape)*2) - np.array([0,101])[:,None,None]
rotcoords = np.dot(np.array([yy,xx]).T, [[np.cos(-38*u.deg), -np.sin(-38*u.deg)],[np.sin(-38*u.deg), np.cos(-38*u.deg)]])
rslt = map_coordinates(nacldata, rotcoords.T)

mn, mid, std = stats.sigma_clipped_stats(nacldata)


pl.clf()
pl.imshow(rslt, origin='lower')
pl.colorbar()

xslice = slice(36,110)
# truncate the disk edges xslice = slice(55,95)
yslice = slice(100,128)
naclcutout = rslt[xslice, yslice]
naclmask = naclcutout > (std*2 + mn)
naclcutout[~naclmask] = np.nan
fig = pl.gcf()
pl.clf()
ax1 = pl.subplot(4,2,1,)
im = ax1.imshow(naclcutout.T, origin='lower',
                extent=[-naclcutout.shape[1]/2 * pixscale * 3600,
                        naclcutout.shape[1]/2 * pixscale * 3600,
                        -naclcutout.shape[0]/2 * pixscale * 3600,
                        naclcutout.shape[0]/2 * pixscale * 3600,
                       ]
                )
ax1.set_aspect(0.15) # hand-picked number, no idea at all where it comes from
#pl.colorbar(mappable=im)
ax2 = pl.subplot(1,2,2)
ax2.plot(np.nanmean(naclcutout, axis=0),
         (np.arange(naclcutout.shape[1]) - naclcutout.shape[1]/2)*pixscale*3600,
         label=linename.replace("v", " v").replace("_"," "))
ax2.set_xlabel("T$_B$ [K]")
ax2.set_ylabel("Offset (\")")
ax2.yaxis.set_label_position('right')
ax2.yaxis.set_ticks_position('right')

print("peaks at {0}".format(scipy.signal.find_peaks(naclcutout.mean(axis=0))[0]))
print("peaks at {0} arcsec".format(scipy.signal.find_peaks(naclcutout.mean(axis=0))[0]*pixscale*3600))
diff = np.diff(scipy.signal.find_peaks(naclcutout.mean(axis=0))[0]*pixscale*3600)
print("height of nacl emission = {0}".format((diff/2*u.arcsec*415*u.pc).to(u.AU, u.dimensionless_angles())))

for ii,contfn in enumerate((b3_hires_cont, b6_hires_cont, b7_hires_cont)):
    conthdu = fits.open(paths.dpath(contfn))[0]

    band = contfn[14:16]
    jtok = radio_beam.Beam.from_fits_header(conthdu.header).jtok(central_freqs[band])

    contcutout,_ = reproject.reproject_interp(conthdu, ww, shape_out=nacldata.shape)
    contcutout *= jtok.value
    reproj_contcutout = map_coordinates(contcutout, rotcoords.T)


    # use sigma clipping to get a mask...
    mn, mid, std = stats.sigma_clipped_stats(contcutout)

    mask = reproj_contcutout > (std*2 + mn)

    reproj_contcutout[~mask] = np.nan

    reproj_cutout = reproj_contcutout[xslice, yslice]

    print("{0}: {1}".format(band, np.nanmean(reproj_cutout[naclmask])))

    ax = pl.subplot(4, 2, 3+ii*2)
    im = ax.imshow(reproj_cutout.T, origin='lower',
                   extent=[-reproj_cutout.shape[1]/2 * pixscale * 3600,
                           reproj_cutout.shape[1]/2 * pixscale * 3600,
                           -reproj_cutout.shape[0]/2 * pixscale * 3600,
                           reproj_cutout.shape[0]/2 * pixscale * 3600,
                          ]
                   )
    ax.set_aspect(0.15) # hand-picked number, no idea at all where it comes from

    ax2.plot(np.nanmean(reproj_cutout, axis=0),
             (np.arange(reproj_cutout.shape[1]) - reproj_cutout.shape[1]/2)*pixscale*3600,
             label="{band} continuum".format(band=band))

ax2.legend(loc='best')
pl.subplot(4,2,1).xaxis.set_ticklabels([])
pl.subplot(4,2,3).xaxis.set_ticklabels([])
pl.subplot(4,2,5).xaxis.set_ticklabels([])
pl.subplot(4,2,7).set_xlabel("Offset (\")")
fig.text(0.02, 0.5, 'Offset (\")', ha='center', va='center', rotation='vertical')

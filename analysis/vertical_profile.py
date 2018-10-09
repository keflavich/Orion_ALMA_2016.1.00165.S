import numpy as np
from spectral_cube import SpectralCube, wcs_utils, tests, Projection
from astropy import units as u
from astropy import wcs
from astropy.io import fits
import paths
import pylab as pl
from scipy.ndimage import map_coordinates
import scipy.signal

sourcename='SourceI'
robust = 0.5
linename='NaClv=2_26-25'
naclfn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.fits').format(linename, sourcename, robust=robust)
nacldata = fits.getdata(naclfn)
ww = wcs.WCS(fits.getheader(naclfn))
pixscale = wcs.utils.proj_plane_pixel_scales(ww)[0]

yy,xx = np.indices(np.array(nacldata.shape)*2) - np.array([0,101])[:,None,None]
rotcoords = np.dot(np.array([yy,xx]).T, [[np.cos(-38*u.deg), -np.sin(-38*u.deg)],[np.sin(-38*u.deg), np.cos(-38*u.deg)]])
rslt = map_coordinates(nacldata, rotcoords.T)

pl.clf(); pl.imshow(rslt, origin='lower'); pl.colorbar()

cutout = rslt[36:110,100:128]
pl.clf()
pl.subplot(1,2,1,)
pl.imshow(cutout.T, origin='lower')
pl.colorbar()
pl.subplot(1,2,2)
pl.plot(cutout.mean(axis=0), np.arange(cutout.shape[1]), )

print("peaks at {0}".format(scipy.signal.find_peaks(cutout.mean(axis=0))[0]))
print("peaks at {0} arcsec".format(scipy.signal.find_peaks(cutout.mean(axis=0))[0]*pixscale*3600))
diff = np.diff(scipy.signal.find_peaks(cutout.mean(axis=0))[0]*pixscale*3600)
print("height of nacl emission = {0}".format((diff/2*u.arcsec*415*u.pc).to(u.AU, u.dimensionless_angles())))

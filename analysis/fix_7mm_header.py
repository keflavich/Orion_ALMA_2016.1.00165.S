"""
Manually shift the Reid+ 2007 7mm data to match ours.  This process invalidates
any sort of absolute astrometry and is eyeballed.
"""
from astropy.io import fits
import paths


fh = fits.open(paths.dpath('srci_reid.7mm.fits'))
fh[0].header['CRPIX1'] = 111.954
fh[0].header['CRPIX2'] = 124.414
fh[0].header['CRVAL1'] = 83.81049
fh[0].header['CRVAL2'] = -5.3751712
fh[0].header['RADESYS'] = 'ICRS'

fh.writeto(paths.dpath('srci_reid.7mm.fixed.fits'), overwrite=True)

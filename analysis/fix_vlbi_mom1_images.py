import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from astropy.table import Table
import regions

import paths

fh = fits.open('/Users/adam/Dropbox/Orion_VLBA/3mm_ORION_CASA_mom1.fits')[0]

# approximate center to be placed on approximate source I center
xcen, ycen = 4182,3793

center_reg = regions.Regions.read(paths.rpath('sourceI_center.reg'))[0]
center = center_reg.center[0]

fits_data = fh.data.squeeze()

ww = wcs.WCS(fh.header).celestial

# fix the central position
ww.wcs.crval[0] = center.ra.deg
ww.wcs.crval[1] = center.dec.deg
ww.wcs.crpix[0] = xcen + 1
ww.wcs.crpix[1] = ycen + 1

good_pixel_indices = np.where(np.isfinite(fits_data))

ra, dec = ww.wcs_pix2world(*good_pixel_indices, 0)*u.deg
vels = fits_data[good_pixel_indices] * u.km/u.s

tbl = Table(data=[ra,dec,vels], names=('RA', 'Dec', 'VLSR'),)
tbl.write(paths.rpath('3mm_maser_velocity_table.fits'), overwrite=True)



# same for 7mm

fh = fits.open('/Users/adam/Dropbox/Orion_VLBA/BG129C_IF1+IF2_mom1.fits')[0]

# approximate center to be placed on approximate source I center
xcen, ycen = 1875,1027

fits_data = fh.data.squeeze()

ww = wcs.WCS(fh.header).celestial

# fix the central position
ww.wcs.crval[0] = center.ra.deg
ww.wcs.crval[1] = center.dec.deg
ww.wcs.crpix[0] = xcen + 1
ww.wcs.crpix[1] = ycen + 1

good_pixel_indices = np.where(np.isfinite(fits_data))

ra, dec = ww.wcs_pix2world(*good_pixel_indices, 0)*u.deg
vels = (fits_data[good_pixel_indices] * u.m/u.s).to(u.km/u.s)

tbl = Table(data=[ra,dec,vels], names=('RA', 'Dec', 'VLSR'),)
tbl.write(paths.rpath('7mm_maser_velocity_table.fits'), overwrite=True)

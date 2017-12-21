import paths
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
from astropy import coordinates, units as u
import radio_beam

coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))

sourceIcont = fits.open(paths.dpath('uid_A001_X88e_X1d3_calibrated_final_cont.pbcor.fits'))

co = Cutout2D(data=sourceIcont[0].data.squeeze(),
              wcs=wcs.WCS(sourceIcont[0].header).celestial, position=coord,
              size=1*u.arcsec)

header = co.wcs.to_header()
beam = radio_beam.Beam.from_fits_header(sourceIcont[0].header)
header.update(beam.to_header_keywords())

hdu = fits.PrimaryHDU(data=co.data, header=header)

hdu.writeto(paths.dpath('OrionSourceI_Band6_QA2_continuum_cutout.fits'),
            overwrite=True)

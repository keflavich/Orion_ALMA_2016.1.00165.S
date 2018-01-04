import paths
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
from astropy import coordinates, units as u
import radio_beam
from radio_beam.beam import NoBeamException

coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))

#sourceIcont = fits.open(paths.dpath('uid_A001_X88e_X1d3_calibrated_final_cont.pbcor.fits'))
for basefn in ('Orion_SourceI_B3_continuum_r-2{suffix}',
               'Orion_SourceI_B6_continuum_r-2_longbaselines{suffix}',
               'Orion_SourceI_B3_continuum_r-2.mask2mJy.clean1mJy{suffix}',
               'Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy{suffix}',
               'Orion_SourceI_B3_continuum_r-2.mask2.5mJy.clean0.5mJy{suffix}',
               'Orion_SourceI_B6_continuum_r0.5{suffix}',
               'Orion_SourceI_B6_continuum_r2{suffix}',
               'uid___A001_X88e_X1df.Orion_BNKL_source_I_sci.spw25_27_29_31.cont{suffix}',
              ):

    for suffix in ('image.tt0.pbcor.fits', 'residual.tt0.fits', 'model.tt0.fits', 'psf.tt0.fits'):

        outfilename = basefn.format(suffix="_SourceIcutout."+suffix)
        fn = basefn.format(suffix="."+suffix)

        sourceIcont = fits.open(paths.dpath(fn))


        co = Cutout2D(data=sourceIcont[0].data.squeeze(),
                      wcs=wcs.WCS(sourceIcont[0].header).celestial, position=coord,
                      size=1*u.arcsec)

        header = co.wcs.to_header()
        try:
            beam = radio_beam.Beam.from_fits_header(sourceIcont[0].header)
            print("Beam is {0} for file {1}".format(beam, fn))
        except NoBeamException:
            print("Assuming beam is the same as previous: {0}".format(beam))
        header.update(beam.to_header_keywords())

        hdu = fits.PrimaryHDU(data=co.data, header=header)

        hdu.writeto(paths.dpath(outfilename),
                    overwrite=True)

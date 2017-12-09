from spectral_cube import SpectralCube
from astropy import units as u

files = [
'/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw0.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw1.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw2.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw3.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw0.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw1.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw2.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw3.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw0.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw1.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw2.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw3.maskedclarkclean10000.image.pbcor.fits',
]

for fn in files:
    cube = SpectralCube.read(fn).mask_out_bad_beams(0.1)

    med = cube.median(axis=0)

    medsub = cube-med

    medsub.write(fn.replace(".image.pbcor.fits", "_medsub.image.pbcor.fits"), overwrite=True)

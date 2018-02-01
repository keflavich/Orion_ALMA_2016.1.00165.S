from spectral_cube import SpectralCube
from astropy import units as u
from astropy import log

files = [
'/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw0.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw1.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw2.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw3.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw0.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw1.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw2.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw3.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw0.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw1.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw2.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw3.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw0.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw1.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw2.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw3.clarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw0.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw1.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw2.maskedclarkclean10000.image.pbcor.fits',
'/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw3.maskedclarkclean10000.image.pbcor.fits',
]

for source in ('SourceI', 'BN'):
    for fn in files:
        try:
            cube = SpectralCube.read(fn.format(source=source))
            if cube.wcs.wcs.radesys.lower() == 'fk5':
                log.exception("{0} has radesys = {1}".format(fn.format(source=source),
                                                             cube.wcs.wcs.radesys))
                continue
            cube = cube.mask_out_bad_beams(0.1)
        except FileNotFoundError:
            log.exception("File {0} not found".format(fn.format(source=source)))
            continue

        log.info("Working on file {0}".format(fn.format(source=source)))

        med = cube.median(axis=0)

        medsub = cube-med

        medsub.write(fn.replace(".image.pbcor.fits",
                                "_medsub.image.pbcor.fits").format(source=source),
                     overwrite=True)

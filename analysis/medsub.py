import os
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import log

import sys
sys.path.append('.')
from source_ids import sources_fmtd

basepath = '/lustre/aginsbur/orion/2016.1.00165.S/imaging'

# files = [
# '/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw0.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw1.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw2.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust0.5.spw3.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw0.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw1.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw2.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B6.robust-2.spw3.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw0.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw1.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw2.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust-2.spw3.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw0.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw1.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw2.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B3.robust0.5.spw3.clarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw0.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw1.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw2.maskedclarkclean10000.image.pbcor.fits',
# '/Volumes/external/orion/Orion{source}_only.B7.robust-2.spw3.maskedclarkclean10000.image.pbcor.fits',
# ]

for source in sources_fmtd:
    for robust in (0.5, 2, -2):
        for suffix, niter in (('maskedclarkclean10000', 10000), ):
            for spw,spws in enumerate([["25","25"], ["27","27"], ["29","29"], ["31","31"], (0,4), (1,5), (2,6), (3,7), (0,), (1,), (2,), (3,)]):
                for band in ("B7.lb", "B6", "B3"):
                    imagename = 'Orion{3}_only.{band}.robust{2}.spw{0}.{1}'.format(spw, suffix, robust, source, band=band)
                    fn = "{basepath}/{imagename}.image.pbcor.fits".format(basepath=basepath, imagename=imagename)
                    doublemedsubfn = fn.replace(".image.pbcor.fits",
                                                "_doublemedsub.image.pbcor.fits").format(source=source)
                    # NOTE: sometimes you want to overwrite stuff!  If so, just do 'if True' here
                    if not os.path.exists(doublemedsubfn):
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

                        medspace = medsub.median(axis=(1,2))
                        doublemedsub = medsub - medspace.quantity[:,None,None]

                        doublemedsub.write(doublemedsubfn,
                                           overwrite=True)


import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
import radio_beam
import paths

files = [
    'OrionSourceI_Unknown_1_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_2_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_3_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_4_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_5_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_8_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_9_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_10_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U229.682_robust0.5maskedclarkclean10000_medsub_K.fits',
]

allbeamlist = [SpectralCube.read(paths.dpath('cubes/'+fn)).beams.common_beam()
               for fn in files]
allbeams = radio_beam.Beams(major=u.Quantity([x.major for x in allbeamlist]),
                            minor=u.Quantity([x.minor for x in allbeamlist]),
                            pa=u.Quantity([x.pa for x in allbeamlist]))
refbeam = allbeams.common_beam()

refcube = SpectralCube.read(paths.dpath('cubes/'+files[0]))
refcube = refcube.convolve_to(refbeam)

stackcube = [refcube]

for cubefn in files[1:]:
    reproj_cube = (SpectralCube.read(paths.dpath('cubes/'+cubefn))
                   .convolve_to(refbeam)
                   .spectral_interpolate(refcube.spectral_axis))

    stackcube.append(reproj_cube)

meancube = np.array([cb.filled_data[:] for cb in stackcube]).mean(axis=0)
avhdu = refcube.hdu
avhdu.data = meancube
avhdu.writeto(paths.dpath('cubes/OrionSourceI_stacked_unknown_lines.fits'))

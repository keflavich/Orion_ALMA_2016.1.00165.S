import numpy as np
import os
import spectral_cube.analysis_utilities
from spectral_cube import SpectralCube
from astropy import units as u
import paths
import pylab as pl

# step 1: create a velocity map

cube = SpectralCube.read(paths.dpath('cubes/OrionSourceI_Unknown_4_robust0.5.maskedclarkclean10000_medsub_K.fits'))
m1 = cube.moment1()
m0 = cube.moment0()
mask = m0.value > 300

vmap = m1
vmap[~mask] = np.nan


# step 2: stack

fullcubepath = '/Volumes/external/orion/'
def fcp(x):
    return os.path.join(fullcubepath, x)

for spw in (0,1,2,3):

    fn = fcp('OrionSourceI_only.B6.robust0.5.spw{0}.maskedclarkclean10000_medsub.image.pbcor.fits').format(spw)
    fullcube = (SpectralCube.read(fn))
    print(fn,fullcube.spectral_extrema)
    fullcube = fullcube.with_spectral_unit(u.km/u.s,
                                           velocity_convention='radio',
                                           rest_value=fullcube.spectral_axis.mean())

    stack = spectral_cube.analysis_utilities.stack_spectra(fullcube, vmap,
                                                           v0=0.0*u.km/u.s)
    fstack = stack.with_spectral_unit(u.GHz)

    fstack.write(paths.dpath('stacked_spectra/OrionSourceI_B6_spw{0}.fits'.format(spw)),
                 overwrite=True)

    pl.clf()
    fstack.quicklook(filename=paths.fpath('stacked_spectra/OrionSourceI_B6_spw{0}.pdf').format(spw))

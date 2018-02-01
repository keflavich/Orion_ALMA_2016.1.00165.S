"""
Tool to produce initial guesses for centroid_planes fitter
"""

import numpy as np
import paths
import regions
from astropy import units as u
from spectral_cube import SpectralCube
from colorize_regions import colorize_region

for fn in (
    'OrionSourceI_29SiOv=0_5-4_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_H2Ov2=1_5(5,0)-6(4,3)_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_H30a_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_HC3N_24-23_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Si34S_13-12_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_SiOv=0_5-4_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_SiOv=1_5-4_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_SiS_12-11_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U217.229_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U218.584_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U229.550_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U229.682_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U229.819_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U230.726_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U230.966_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_U232.634_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_10_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_11_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_12_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_14_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_15_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_16_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_1_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_2_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_3_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_4_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_5_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_6_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_7_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_8_robust0.5maskedclarkclean10000_medsub_K.fits',
    'OrionSourceI_Unknown_9_robust0.5maskedclarkclean10000_medsub_K.fits',
):
    cube = SpectralCube.read(paths.dpath('cubes/{0}'.format(fn)))

    std = cube.mad_std(axis=(1,2))
    mx = cube.max(axis=(1,2))

    mask = (mx > std*5)
    print("Found {0} entries for {1}".format(mask.sum(), fn))
    if mask.sum() < 2:
        continue

    rname = fn.replace("robust0.5maskedclarkclean10000_medsub_K.fits", "initial_guesses.reg")
    rpath = paths.rpath('initial_guesses_centroiding/'+rname)
    with open(rpath, 'w') as fh:
        fh.write('image\n')
        for ind in np.where(mask)[0]:
            ymx, xmx = np.unravel_index(cube[int(ind),:,:].argmax(), cube.shape[1:])
            vel = cube.spectral_axis[int(ind)]

            fh.write("point({0}, {1}) # point=x text={{{2}}}\n"
                     .format(xmx+1, ymx+1, np.round(vel.to(u.km/u.s).value, 2)))

    regs = colorize_region(rpath)
    regions.io.ds9.write_ds9(regs, filename=rpath, coordsys='image')

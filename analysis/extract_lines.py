import paths
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from spectral_cube import SpectralCube, wcs_utils, tests

ftemplate = '/Volumes/external/orion/Orion{1}_only.B6.robust0.5.spw{0}.maskedclarkclean10000_medsub.image.pbcor.fits'

for sourcename in ('SourceI', 'BN'):
    for linename, linefreq, vrange in [('H2Ov2=1_5(5,0)-6(4,3)', 232.6867*u.GHz, [-60, 70]),
                                       ('H30a', 231.900928*u.GHz, [-120, 130]),
                                       ('Unknown_1', 230.321535*u.GHz, [-60,70]),
                                       ('SiS_12-11', 217.817644*u.GHz, [-60,70]),
                                      ]:
        for spw in (0,1,2,3):
            cube = SpectralCube.read(ftemplate.format(spw, sourcename))

            fmin, fmax = cube.spectral_extrema

            if linefreq > fmin and linefreq < fmax:

                print("Writing {0} in spw {1}".format(linename, spw))

                scube = (cube.with_spectral_unit(u.km/u.s,
                                                 velocity_convention='radio',
                                                 rest_value=linefreq)
                         .spectral_slab(-60*u.km/u.s, 70*u.km/u.s))

                cubeK = scube.to(u.K)
                cubeK.write(paths.dpath('cubes/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K.fits').format(linename, sourcename),
                            overwrite=True)

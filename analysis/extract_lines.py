import paths
import os
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from spectral_cube import SpectralCube, wcs_utils, tests
import pylab as pl
from lines import disk_lines

ftemplate = '/Volumes/external/orion/Orion{1}_only.B6.robust0.5.spw{0}.maskedclarkclean10000_medsub.image.pbcor.fits'

conthdu = fits.open(paths.dpath('OrionSourceI_Band6_QA2_continuum_cutout.fits'))

for sourcename in ('SourceI', 'BN'):
    for linename, linefreq in disk_lines.items():
        if 'H30a' in linename:
            vrange = [-120,130]
        else:
            vrange = [-30,40]

        for spw in (0,1,2,3):
            filename = ftemplate.format(spw, sourcename)
            if os.path.exists(filename):
                cube = SpectralCube.read(filename)
            else:
                log.exception("File {0} does not exist".format(filename))
                continue

            fmin, fmax = cube.spectral_extrema

            if linefreq > fmin and linefreq < fmax:

                print("Writing {0} in spw {1}".format(linename, spw))

                scube = (cube.with_spectral_unit(u.km/u.s,
                                                 velocity_convention='radio',
                                                 rest_value=linefreq)
                         .spectral_slab(vrange[0]*u.km/u.s,
                                        vrange[1]*u.km/u.s))

                cubeK = scube.to(u.K)
                cubeK.write(paths.dpath('cubes/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K.fits').format(linename, sourcename),
                            overwrite=True)

                m0 = cubeK.moment0(axis=0)
                mx = cubeK.max(axis=0)

                mx.write(paths.dpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_peak.fits').format(linename,
                                                                                                                       sourcename),
                         overwrite=True)
                m0.write(paths.dpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_moment0.fits').format(linename,
                                                                                                                          sourcename),
                         overwrite=True)

                mx.quicklook(filename=paths.fpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_peak.png')
                             .format(linename, sourcename), aplpy_kwargs={'figure':pl.figure(1)})
                mx.FITSFigure.show_contour(conthdu, levels=[0.001, 0.005, 0.01,
                                                            0.02, 0.03, 0.04,
                                                            0.05],
                                           colors=['r']*10)
                mx.FITSFigure.colorbar.set_axis_label_text("$T_B$ [K]")
                mx.FITSFigure.save(paths.fpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_peak.png')
                                   .format(linename, sourcename))

                pl.figure(1).clf()

                m0.quicklook(filename=paths.fpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_moment0.png')
                             .format(linename, sourcename), aplpy_kwargs={'figure':pl.figure(1)})

                m0.FITSFigure.show_contour(conthdu, levels=[0.001, 0.005, 0.01,
                                                            0.02, 0.03, 0.04,
                                                            0.05],
                                           colors=['r']*10)
                m0.FITSFigure.colorbar.set_axis_label_text("$\int T_B \mathrm{d}v$ [K km s$^{-1}$]")
                m0.FITSFigure.save(paths.fpath('moments/Orion{1}_{0}_robust0.5.maskedclarkclean10000_medsub_K_moment0.png')
                                   .format(linename, sourcename))

                pl.figure(1).clf()

import paths
import os
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from spectral_cube import SpectralCube, wcs_utils, tests
import pylab as pl

ftemplate = '/Volumes/external/orion/Orion{1}_only.{2}.robust0.5.spw{0}.{suffix}_medsub.image.pbcor.fits'

conthdu = fits.open(paths.dpath('OrionSourceI_Band6_QA2_continuum_cutout.fits'))

for band,suffix in (('B3', 'clarkclean10000'),
                    ('B6', 'maskedclarkclean10000')):
    for sourcename in ('SourceI', 'BN'):

        for spw in (0,1,2,3):
            filename = ftemplate.format(spw, sourcename, band, suffix=suffix)
            if os.path.exists(filename):
                cube = SpectralCube.read(filename)
            else:
                log.exception("File {0} does not exist".format(filename))
                continue


            outf_template = os.path.basename(filename).replace(".image.pbcor","_{0}")

            print("Writing in spw {0}".format(spw))

            cube.allow_huge_operations = True
            cubeK = cube.to(u.K)
            avgspec = cubeK.mean(axis=(1,2))
            stdspec = cubeK.std(axis=(1,2))
            maxspec = cubeK.max(axis=(1,2))
            minspec = cubeK.min(axis=(1,2))

            avgspec.write(paths.dpath('spectra/{0}'.format(outf_template.format('avg'))), overwrite=True)
            stdspec.write(paths.dpath('spectra/{0}'.format(outf_template.format('std'))), overwrite=True)
            minspec.write(paths.dpath('spectra/{0}'.format(outf_template.format('min'))), overwrite=True)
            maxspec.write(paths.dpath('spectra/{0}'.format(outf_template.format('max'))), overwrite=True)

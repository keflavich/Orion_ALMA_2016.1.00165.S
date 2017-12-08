import numpy as np
import os
import glob
import paths
from astropy import table
from astropy import units as u
from astropy import constants
from astropy import coordinates
from astropy import wcs
from astropy import log
from astroquery.splatalogue import Splatalogue
from spectral_cube import SpectralCube
import pyregion
import pyspeckit
import imp; import lines; imp.reload(lines)
from lines import disk_lines

# use outflow_meta b/c higher precision than ds9 reg

import pylab as pl

pl.close(1)

diskycoorddict = {}
source = "sourceI"
coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))
extraction_region1, extraction_region2 = pyregion.open(paths.rpath('spectral_extraction_regions.reg'))

linedata_fits = {}

for name, cutoutname, source, vrange, vcen in (
    ('sourceI', 'sourceI', coord, (-30,40), 6.5),
   ):

    for spw in (0,1,2,3):

        vcen = u.Quantity(vcen, u.km/u.s)

        #fn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout.fits'
        fn = '/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw{0}.maskedclarkclean10000.image.pbcor.fits'.format(spw)

        #cube = (SpectralCube.read(fn)[:,515:721,550:714].mask_out_bad_beams(5))
        cube = (SpectralCube.read(fn).mask_out_bad_beams(5))
        # cube.allow_huge_operations=True
        cube.beam_threshold = 5000
        log.info("Calculating 25th percentile")
        med = cube.percentile(25,axis=0)
        medsub = cube - med

        for linename, linefreq in disk_lines.items():

            subcube = (medsub.with_spectral_unit(u.km/u.s,
                                                 velocity_convention='radio',
                                                 rest_value=linefreq)
                       .spectral_slab(-50*u.km/u.s, 60*u.km/u.s)
                       .with_spectral_unit(u.GHz))

            if subcube.shape[0] < 5:
                log.warn("Skipping line {0} in {1} because it's empty".format(linename, fn))
                continue

            specblue = subcube.subcube_from_ds9region(pyregion.ShapeList([extraction_region1])).mean(axis=(1,2))
            specred = subcube.subcube_from_ds9region(pyregion.ShapeList([extraction_region2])).mean(axis=(1,2))

            spred = pyspeckit.Spectrum.from_hdu(specred.hdu.copy())
            spblue = pyspeckit.Spectrum.from_hdu(specblue.hdu.copy())

            assert spred.xarr.unit == u.GHz

            vred = vcen+12*u.km/u.s
            vblue = vcen-12*u.km/u.s

            spred.specfit(guesses=[spred.data.max(), (1-vred/constants.c)*linefreq.value, 2/3e5*linefreq.value])
            spblue.specfit(guesses=[spblue.data.max(), (1-vblue/constants.c)*linefreq.value, 2/3e5*linefreq.value])

            avfrq = (spred.specfit.parinfo[1].value + spblue.specfit.parinfo[1].value)/2. * u.GHz
            corrected_avfrq = avfrq * (1+vcen/constants.c)
            vshift = (corrected_avfrq-linefreq)/linefreq * constants.c.to(u.km/u.s)
            print("For line {0}, the average frequency is {1} (orig was {2})".format(linename, corrected_avfrq, linefreq))
            print("Shift is {0} - difference is {1}".format(vshift, vshift))
    
            v = vcen
            dv = 2*u.km/u.s
            frq = avfrq
            result = Splatalogue.query_lines(frq + (v-dv)/constants.c*frq, frq + (v+dv)/constants.c*frq)
            ref=np.array(result['Freq-GHz'])*u.GHz
            result.add_column(table.Column(name='velocity', data=-((frq-ref)/(ref) * constants.c).to(u.km/u.s)))

            linedata_fits[linename] = {'avfrq': corrected_avfrq,
                                       'redfit': spred.specfit.parinfo,
                                       'bluefit': spblue.specfit.parinfo,
                                       'splat': result['Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','velocity', 'E_U (K)'],
                                       'vshift': vshift,
                                       'linename': linename,
                                       'linefreq': linefreq,
                                      }

            #pl.figure(1)
            #spblue.plotter(axis=pl.subplot(2,1,1))
            #spblue.specfit.plot_fit()
            #spred.plotter(axis=pl.subplot(2,1,2))
            #spred.specfit.plot_fit()

            #input("")

for entry in linedata_fits:
    print("Line {linename} is shifted {vshift} and should be {avfrq}".format(**linedata_fits[entry]))

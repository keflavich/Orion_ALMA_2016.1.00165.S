import numpy as np
import pvextractor
import os
import glob
import paths
from astropy import units as u
from astropy import constants
from astropy import coordinates
from astropy import wcs
from astropy import log
from spectral_cube import SpectralCube
import pyregion
import imp; import lines; imp.reload(lines)
from lines import disk_lines
import show_pv; imp.reload(show_pv)
import re

robustnumre = re.compile('robust(0.5|-2|2)')

# use outflow_meta b/c higher precision than ds9 reg
from line_point_offset import offset_to_point

import pylab as pl

pl.close(1)

diskycoorddict = {}
source = "sourceI"
coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))
diskycoord_list = pyregion.open(paths.rpath("{0}_disk_pvextract.reg"
                                            .format(source)))[0].coord_list
diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoord_list[jj],
                                                     diskycoord_list[jj+1])
                                    for jj in range(0,
                                                    len(diskycoord_list),
                                                    2)], unit=(u.deg,
                                                               u.deg),
                                   frame='fk5')
diskycoorddict[source] = diskycoords

for width in (0.1, 0.01):
    for name, cutoutname, source, vrange, vcen in (
        ('sourceI', 'sourceI', coord, (-30,40), 6.0),
       ):

        for fnt in ('/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    '/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw{0}.clarkclean10000.image.pbcor.fits',
                    '/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw{0}.clarkclean10000.image.pbcor.fits',
                   ):

            for spw in (0,1,2,3):

                fn = fnt.format(spw)

                vcen = u.Quantity(vcen, u.km/u.s)

                #fn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout.fits'

                #cube = (SpectralCube.read(fn)[:,515:721,550:714].mask_out_bad_beams(5))
                cube = (SpectralCube.read(fn).mask_out_bad_beams(5))
                # cube.allow_huge_operations=True
                cube.beam_threshold = 5000
                log.info("Calculating 25th percentile")
                med = cube.percentile(25,axis=0)
                medsub = cube - med

                # width = 0.05 arcsec encompasses the disk; however, most
                # of the actual line emission comes from just above/below...
                #extraction_path = pvextractor.Path(diskycoords, width=0.05*u.arcsec)
                extraction_path = pvextractor.Path(diskycoords, width=width*u.arcsec)
                log.info("Beginning extraction of path with width {0}".format(extraction_path.width))
                extracted = pvextractor.extract_pv_slice(medsub, extraction_path)
                outfn = paths.dpath(os.path.join('pv',
                                                 os.path.split(fn.replace(".image.pbcor.fits",
                                                                          "_medsub_diskpv_{0}.fits".format(width)))[-1]))
                log.info("Writing to {0}".format(outfn))
                extracted.writeto(outfn, overwrite=True)

                for linename, linefreq in disk_lines.items():

                    diskycoords = diskycoorddict[name]

                    subcube = (medsub.with_spectral_unit(u.km/u.s,
                                                         velocity_convention='radio',
                                                         rest_value=linefreq)
                               .spectral_slab(-50*u.km/u.s, 60*u.km/u.s))

                    if subcube.shape[0] < 5:
                        log.warn("Skipping line {0} in {1} because it's empty".format(linename, fn))
                        continue

                    if 'robust' in fn:
                        robustnum = robustnumre.search(fn).groups()[0]
                        basename = "{0}_{1}_robust{2}_diskpv_{3}.fits".format(name,
                                                                              linename,
                                                                              robustnum,
                                                                              width)
                    else:
                        basename = "{0}_{1}_diskpv_{2}.fits".format(name, linename, width)
                    outfn = paths.dpath(os.path.join("pv/", basename))

                    extraction_path = pvextractor.Path(diskycoords, width*u.arcsec)
                    log.info("Beginning extraction of path with width {0}".format(extraction_path.width))
                    extracted = pvextractor.extract_pv_slice(subcube, extraction_path)
                    log.info("Writing to {0}".format(outfn))
                    extracted.writeto(outfn, overwrite=True)

                    origin = offset_to_point(source.ra.deg,
                                             source.dec.deg,
                                             extraction_path)*u.deg

                    ww = wcs.WCS(extracted.header)
                    ww.wcs.cdelt[1] /= 1000.0
                    ww.wcs.crval[1] /= 1000.0
                    ww.wcs.cunit[1] = u.km/u.s
                    ww.wcs.cdelt[0] *= 3600
                    ww.wcs.cunit[0] = u.arcsec
                    ww.wcs.crval[0] = -origin.to(u.arcsec).value

                    # #!@#$!@#$@!%@#${^(@#$)%#$(
                    ww.wcs.set()

                    if ww.wcs.cunit[1] == 'm/s':
                        scalefactor = 1000.0
                    else:
                        scalefactor = 1.0

                    good_limits = (np.array((np.argmax(np.isfinite(extracted.data.max(axis=0))),
                                             extracted.data.shape[1] -
                                             np.argmax(np.isfinite(extracted.data.max(axis=0)[::-1])) - 1
                                            ))
                                   )
                    leftmost_position = ww.wcs_pix2world(good_limits[0],
                                                         vrange[0]*scalefactor,
                                                         0)[0]*u.arcsec
                    rightmost_position = ww.wcs_pix2world(good_limits[1],
                                                          vrange[0]*scalefactor,
                                                          0)[0]*u.arcsec
                    assert rightmost_position > 0
                    maxdist = ((rightmost_position)*415*u.pc).to(u.au, u.dimensionless_angles())
                    assert maxdist > 0

                    plotted_region = ww.wcs_world2pix([0,0],
                                                      np.array(vrange)*scalefactor,
                                                      0)
                    plotted_slice = (slice(int(np.min(plotted_region[1])), int(np.max(plotted_region[1]))),
                                     slice(None,None),
                                    )
                    vmin,vmax = (np.nanmin(extracted.data[plotted_slice]),
                                 np.nanmax(extracted.data[plotted_slice]))
                    vmin = -0.0025
                    if vmax < 0.5:
                        vmax = np.min([0.02, vmax])
                    if 'H2O' in linename and width > 0.05:
                        vmax = 0.1

                    fig,ax = show_pv.show_pv(extracted.data, ww,
                                             origin, vrange=vrange, vcen=vcen,
                                             imvmin=vmin, imvmax=vmax)


                    ax.set_xlim(good_limits)

                    fig.savefig(paths.fpath('pv/{0}/'.format(name, linename) +
                                            basename.replace(".fits",".png")),
                                dpi=200,
                                bbox_inches='tight')


                    # override that previous junk since we went through the effort of calculating it
                    ax.set_xlim(good_limits)

                    show_pv.show_keplercurves(ax, origin, maxdist, vcen,
                                              masses=[19],
                                              linestyles='-',
                                             )

                    fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name, linename) +
                                            basename.replace(".fits",".png")),
                                dpi=200,
                                bbox_inches='tight')

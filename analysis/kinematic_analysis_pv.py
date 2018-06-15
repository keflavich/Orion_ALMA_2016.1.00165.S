import numpy as np
import pvextractor
import os
import glob
import paths
import json
from astropy.io import fits
from astropy import units as u
from astropy import constants
from astropy import stats
from astropy import coordinates
from astropy import wcs
from astropy import log
from astropy.table import Table
from spectral_cube import SpectralCube
import regions
import imp; import lines; imp.reload(lines)
from lines import disk_lines
import show_pv; imp.reload(show_pv)
import re
from mpl_plot_templates import adaptive_param_plot
from constants import vcen as assumed_vcen
import time

robustnumre = re.compile('robust(0.5|-2|2)')
bandre = re.compile("\.B([367])\.")

# use outflow_meta b/c higher precision than ds9 reg
from line_point_offset import offset_to_point

import pylab as pl

pl.ioff()

pl.close(1)

redo = False

errtbl = {}

t0 = time.time()

# just do water (goes faster)
#disk_lines = {x:y for x,y in disk_lines.items() if 'H2O' in x}
#disk_lines = {x:y for x,y in disk_lines.items() if 'H2O' in x or 'SiO' in x}
#disk_lines = {x:y for x,y in disk_lines.items() if 'Unknown_4' in x or 'Unknown_5' in x}
#disk_lines = {x:y for x,y in disk_lines.items() if 'U230.966' in x or 'U230.72' in x}
#redo = True

diskycoorddict = {}
source = "sourceI"
coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))

# this is the fitted disk center value from the disk continuum modeling
# This position *must* be in the same frame as diskycoordlist below
coord = coordinates.SkyCoord(83.81048816210084, -5.3751716623649575,
                             frame='icrs',
                             unit=(u.hour, u.deg))

#diskycoord_list = pyregion.open(paths.rpath("{0}_disk_pvextract.reg"
#                                            .format(source)))[0].coord_list
#diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoord_list[jj],
#                                                     diskycoord_list[jj+1])
#                                    for jj in range(0,
#                                                    len(diskycoord_list),
#                                                    2)], unit=(u.deg,
#                                                               u.deg),
#                                   frame='icrs')
#diskycoorddict[source] = diskycoords
diskycoord_list = regions.read_ds9(paths.rpath("{0}_disk_pvextract.reg"
                                               .format(source)))
diskycoorddict[source] = coordinates.SkyCoord([diskycoord_list[0].start,
                                               diskycoord_list[0].end])


for width in (0.01, 0.05, 0.1, 0.2, 0.3, 0.4):
    for name, cutoutname, source, vrange, vcen in (
        ('sourceI', 'sourceI', coord, (-30,40), assumed_vcen),
       ):

        diskycoords = diskycoorddict[name]

        for fnt in (
                    #'/Volumes/external/orion/OrionSourceI_only.B6.robust0.5.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    #'/Volumes/external/orion/OrionSourceI_only.B6.robust-2.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    #'/Volumes/external/orion/OrionSourceI_only.B6.robust-2.longbaselines.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    #'/Volumes/external/orion/OrionSourceI_only.B3.robust-2.spw{0}.clarkclean10000.image.pbcor.fits',
                    #'/Volumes/external/orion/OrionSourceI_only.B3.robust0.5.spw{0}.clarkclean10000.image.pbcor.fits',
                    '/Volumes/external/orion/OrionSourceI_only.B7.lb.robust0.5.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    '/Volumes/external/orion/OrionSourceI_only.B7.lb.robust-2.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                    '/Volumes/external/orion/OrionSourceI_only.B7.lb.robust2.spw{0}.maskedclarkclean10000.image.pbcor.fits',
                   ):

            for spw in (0,1,2,3):


                if 'longbaselines' in name:
                    name = name+"_longbaselines"

                fn = fnt.format(spw)

                vcen = u.Quantity(vcen, u.km/u.s)

                #fn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout.fits'

                medsubfn = fn.replace(".image.pbcor.fits",
                                      "_medsub.image.pbcor.fits")

                if os.path.exists(medsubfn):
                    medsub = SpectralCube.read(medsubfn)
                    medsub.beam_threshold = 5000
                    if not medsub.wcs.wcs.radesys.lower() == 'icrs':
                        log.exception("Skipping {0} because of a bad coordinate system.".format(medsubfn))
                        continue

                else:
                    #cube = (SpectralCube.read(fn)[:,515:721,550:714].mask_out_bad_beams(5))
                    cube = (SpectralCube.read(fn).mask_out_bad_beams(5))
                    if not cube.wcs.wcs.radesys.lower() == 'icrs':
                        log.exception("Skipping {0} because of a bad coordinate system.".format(fn))
                        continue
                    # cube.allow_huge_operations=True
                    cube.beam_threshold = 5000
                    log.info("Calculating 25th percentile")
                    med = cube.percentile(25,axis=0)
                    medsub = cube - med
                    medsub.write(medsubfn)


                # create full-scale PV diagram (all freqs)
                outfn = paths.dpath(os.path.join('pv',
                                                 os.path.split(fn.replace(".image.pbcor.fits",
                                                                          "_medsub_diskpv_{0}.fits".format(width)))[-1]))
                if not os.path.exists(outfn):
                    # width = 0.05 arcsec encompasses the disk; however, most
                    # of the actual line emission comes from just above/below...
                    #extraction_path = pvextractor.Path(diskycoords, width=0.05*u.arcsec)
                    extraction_path = pvextractor.Path(diskycoords, width=width*u.arcsec)
                    log.info("Beginning extraction of path with width {0} for {1}".format(extraction_path.width, outfn))
                    extracted = pvextractor.extract_pv_slice(medsub, extraction_path)
                    log.info("Writing to {0}".format(outfn))
                    extracted.writeto(outfn, overwrite=True)


                for linename, linefreq in disk_lines.items():
                    print(linename, width, fnt, spw, time.time() - t0)


                    pl.close('all')

                    band = 'B'+bandre.search(fnt).groups()[0]

                    if 'robust' in fn:
                        robustnum = robustnumre.search(fn).groups()[0]
                        basename = ("{0}_{1}_{4}_robust{2}_diskpv_{3}.fits"
                                    .format(name, linename, robustnum, width, band)
                                   )
                        robustnum = float(robustnum)
                    else:
                        robustnum = 999
                        basename = ("{0}_{1}_{3}_diskpv_{2}.fits"
                                    .format(name, linename, width, band))
                    outfn = paths.dpath(os.path.join("pv/", basename))

                    extraction_path = pvextractor.Path(diskycoords, width*u.arcsec)

                    if os.path.exists(outfn) and not redo:
                        extracted = fits.open(outfn)[0]
                        if 'RADESYS' in extracted.header and not extracted.header['RADESYS'].lower == 'icrs':
                            log.exception("Skipping line {0} in {1} because of a bad coordinate system.".format(linename, fn))
                            continue
                    else:
                        subcube = (medsub.with_spectral_unit(u.km/u.s,
                                                             velocity_convention='radio',
                                                             rest_value=linefreq)
                                   .spectral_slab(-50*u.km/u.s, 60*u.km/u.s))

                        if subcube.shape[0] < 5:
                            log.warn("Skipping line {0} in {1} because it's empty".format(linename, fn))
                            continue


                        log.info("Beginning extraction of path with width {0} for {1}".format(extraction_path.width, outfn))
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

                    ww.wcs.crval[0] = 0
                    ww.wcs.crpix[0] = extracted.data.shape[1]/2+1
                    origin = 0*u.arcsec

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
                    if vmax < 0.5 and 'B7' not in fnt:
                        vmax = np.min([0.02, vmax])
                    if 'H2O' in linename and width > 0.05 and robustnum > -2:
                        vmax = 0.1

                    fig,ax,cb,con = show_pv.show_pv(extracted.data, ww, origin,
                                                    vrange=vrange, vcen=vcen,
                                                    imvmin=vmin, imvmax=vmax,
                                                    contour='b'
                                                   )
                    errtbl[(name, linename, width, robustnum)] = stats.mad_std(extracted.data, ignore_nan=True)


                    ax.set_xlim(good_limits)

                    cb.set_label("$S_\\nu$ [Jy beam$^{-1}$]")

                    bb = ax.bbox._bbox
                    cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])

                    fig.savefig(paths.fpath('pv/{0}/'.format(name, linename) +
                                            basename.replace(".fits","_withcontour.pdf")),
                                dpi=200,
                                bbox_inches='tight')

                    con.set_alpha(0)

                    fig.savefig(paths.fpath('pv/{0}/'.format(name, linename) +
                                            basename.replace(".fits",".pdf")),
                                dpi=200,
                                bbox_inches='tight')

                    # override that previous junk since we went through the effort of calculating it
                    ax.set_xlim(good_limits)

                    show_pv.show_keplercurves(ax, origin, maxdist, vcen,
                                              masses=[15, ],
                                              radii={15: ([30, 80], ['m', 'm'])},
                                              linestyles='-',
                                              colors=['r'],
                                             )

                    bb = ax.bbox._bbox
                    cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])


                    fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name, linename) +
                                            basename.replace(".fits",".pdf")),
                                dpi=200,
                                bbox_inches='tight')

                    bb = ax.bbox._bbox
                    cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])

                    con.set_alpha(1)
                    for coll in con.collections:
                        coll.set_color('w')
                    fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name, linename) +
                                            basename.replace(".fits","_withcontour.pdf")),
                                dpi=200,
                                bbox_inches='tight')




diskycoords = diskycoorddict['sourceI']
extraction_path = pvextractor.Path(diskycoords)

masers_3mm = Table.read(paths.rpath('3mm_maser_velocity_table.fits'))
masers_7mm = Table.read(paths.rpath('7mm_maser_velocity_table.fits'))


# this is the center position used to reference the maser positions
maser_center_reg = regions.read_ds9(paths.rpath('sourceI_center.reg'))[0]
maser_center = maser_center_reg.center

xpoints_3mm = u.Quantity(list(map(lambda x,y:
                                  offset_to_point(x,y,extraction_path),
                                  masers_3mm['RA'], masers_3mm['Dec'])), u.deg)
xpoints_7mm = u.Quantity(list(map(lambda x,y:
                                  offset_to_point(x,y,extraction_path),
                                  masers_7mm['RA'], masers_7mm['Dec'])), u.deg)
maser_center_3mm = (xpoints_3mm.max() + xpoints_3mm.min())/2.
maser_center_7mm = xpoints_7mm.mean()


for owidth,iwidth in ((0.1,0.01), (0.2,0.1), (0.3,0.2), (0.2,0.05), (0.4,0.3)):
    for name, cutoutname, source, vrange, vcen in (
        ('sourceI', 'sourceI', coord, (-30,40), assumed_vcen),
       ):
        for robustnum in (0.5, -2):
            for linename, linefreq in disk_lines.items():

                inner_fn = "{0}_{1}_B6_robust{2}_diskpv_{3}.fits".format(name,
                                                                         linename,
                                                                         robustnum,
                                                                         iwidth)
                outer_fn = "{0}_{1}_B6_robust{2}_diskpv_{3}.fits".format(name,
                                                                         linename,
                                                                         robustnum,
                                                                         owidth)
                try:
                    outerfh = fits.open(paths.dpath(os.path.join("pv/", outer_fn)))
                    innerfh = fits.open(paths.dpath(os.path.join("pv/", inner_fn)))
                except:
                    print("Skipping {0}".format(outer_fn))
                    continue

                pixscale = innerfh[0].header['CDELT1']*u.deg

                outerarea = pixscale * u.Quantity(owidth, u.arcsec)
                innerarea = pixscale * u.Quantity(iwidth, u.arcsec)

                diff = ((outerfh[0].data * outerarea -
                         innerfh[0].data * innerarea) /
                        (outerarea-innerarea)).decompose().value

                outfn = "{0}_{1}_robust{2}_diskpv_{3}-{4}.fits".format(name,
                                                                       linename,
                                                                       robustnum,
                                                                       owidth,
                                                                       iwidth)

                outerfh[0].data = diff
                outerfh.writeto(outfn, overwrite=True)

                ww = wcs.WCS(outerfh[0].header)
                ww.wcs.cdelt[1] /= 1000.0
                ww.wcs.crval[1] /= 1000.0
                ww.wcs.cunit[1] = u.km/u.s
                ww.wcs.cdelt[0] *= 3600
                ww.wcs.cunit[0] = u.arcsec

                # #!@#$!@#$@!%@#${^(@#$)%#$(
                ww.wcs.set()

                if ww.wcs.cunit[1] == 'm/s':
                    scalefactor = 1000.0
                else:
                    scalefactor = 1.0

                #origin_ = ww.sub([1]).all_pix2world([outerfh[0].data.shape[1]/2], 0)[0][0]*u.arcsec
                ww.wcs.crval[0] = 0
                ww.wcs.crpix[0] = outerfh[0].data.shape[1]/2+1
                origin = 0*u.arcsec

                vmin,vmax = (np.nanmin(diff),
                             np.nanmax(diff))
                vmin = -0.0025
                if vmax < 0.5 and 'B7' not in outfn:
                    vmax = np.min([0.02, vmax])
                if 'H2O' in linename and owidth > 0.05 and robustnum > -2:
                    vmax = 0.1

                fig,ax,cb,con = show_pv.show_pv(diff, ww, origin,
                                                vrange=vrange,
                                                vcen=u.Quantity(vcen,u.km/u.s),
                                                imvmin=vmin, imvmax=vmax,
                                                contour='w',
                                               )
                errtbl[(name, linename, (iwidth, owidth), robustnum)] = stats.mad_std(diff, ignore_nan=True)
                cb.set_label("$S_\\nu$ [Jy beam$^{-1}$]")
                bb = ax.bbox._bbox
                cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])
                fig.savefig(paths.fpath('pv/{0}/'.format(name) +
                                        outfn.replace(".fits","_withcontour.pdf")),
                            dpi=200,
                            bbox_inches='tight')

                con.set_alpha(0)

                kc = show_pv.show_keplercurves(ax, origin, 150*u.au, u.Quantity(vcen,u.km/u.s),
                                               masses=[15,],
                                               linestyles='-',
                                               colors=['r'],
                                              )

                bb = ax.bbox._bbox
                cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])

                fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name) +
                                        outfn.replace(".fits",".pdf")),
                            dpi=200,
                            bbox_inches='tight')

                con.set_alpha(1)

                fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name) +
                                        outfn.replace(".fits","_withcontour.pdf")),
                            dpi=200,
                            bbox_inches='tight')


                if not ('SiO' in linename or 'H2O' in linename):
                    continue

                for line in kc:
                    line.set_visible(False)


                trans = ax.get_transform('world')
                m3m = ax.plot((xpoints_3mm-maser_center_3mm).to(u.arcsec),
                              u.Quantity(masers_3mm['VLSR'], u.m/u.s), 'r,',
                              transform=trans)
                fig.savefig(paths.fpath('pv/{0}/vlba_3mm_maseroverlay_'.format(name) +
                                        outfn.replace(".fits",".pdf")),
                            dpi=200,
                            bbox_inches='tight')

                for line in m3m:
                    line.set_visible(False)

                trans = ax.get_transform('world')
                #m7m = ax.contour((xpoints_7mm-maser_center_7mm).to(u.arcsec),
                #                 u.Quantity(masers_7mm['VLSR'], u.m/u.s),
                #                 levels=[5,10,50],
                #                 transform=trans)
                m7m = adaptive_param_plot((xpoints_7mm-maser_center_7mm).to(u.arcsec).value,
                                          u.Quantity(masers_7mm['VLSR'], u.m/u.s).value,
                                          axis=ax,
                                          bins=100,
                                          marker=',',
                                          marker_color='r',
                                          transform=trans)


                fig.savefig(paths.fpath('pv/{0}/vlba_7mm_maseroverlay_'.format(name) +
                                        outfn.replace(".fits",".pdf")),
                            dpi=200,
                            bbox_inches='tight')

errtbl2 = {key[1]: (key[2], key[3], value) for key, value in errtbl.items()}
with open('error_estimate_table_pvs.json', 'w') as fh:
    json.dump(errtbl2, fh)

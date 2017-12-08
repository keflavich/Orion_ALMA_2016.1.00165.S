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

        # width = 0.05 arcsec encompasses the disk; however, most
        # of the actual line emission comes from just above/below...
        #extraction_path = pvextractor.Path(diskycoords, width=0.05*u.arcsec)
        extraction_path = pvextractor.Path(diskycoords, width=0.1*u.arcsec)
        log.info("Beginning extraction")
        extracted = pvextractor.extract_pv_slice(medsub, extraction_path)
        outfn = paths.dpath(os.path.join('pv',
                                         os.path.split(fn.replace(".image.pbcor.fits",
                                                                  "_medsub_diskpv.fits"))[-1]))
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


            basename = "{0}_{1}_diskpv.fits".format(name, linename)
            outfn = paths.dpath(os.path.join("pv/", basename))

            extraction_path = pvextractor.Path(diskycoords, 0.05*u.arcsec)
            log.info("Beginning extraction")
            extracted = pvextractor.extract_pv_slice(subcube, extraction_path)
            log.info("Writing to {0}".format(outfn))
            extracted.writeto(outfn, overwrite=True)

            ww = wcs.WCS(extracted.header)
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

            plotted_region = ww.wcs_world2pix([0,0],
                                              np.array(vrange)*scalefactor,
                                              0)
            plotted_slice = (slice(int(np.min(plotted_region[1])), int(np.max(plotted_region[1]))),
                             slice(None,None),
                            )
            if np.any(np.array(extracted.data[plotted_slice].shape) == 0):
                log.warn("Skipping {0} because it's empty".format(fn))
                continue


            fig = pl.figure(1, figsize=(12,8))
            fig.clf()
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8],projection=ww)
            assert ww.wcs.cunit[1] == 'm/s' # this is BAD BAD BAD but necessary

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


            vmin,vmax = (np.nanmin(extracted.data[plotted_slice]),
                         np.nanmax(extracted.data[plotted_slice]))
            vmin = -0.0025
            vmax = np.min([0.02, vmax])
            im = ax.imshow(extracted.data, cmap='gray_r',
                           vmin=vmin, vmax=vmax*1.1,
                           interpolation='none')
            ax.set_xlabel("Offset [\"]")
            ax.set_ylabel("$V_{LSR}$ [km/s]")


            trans = ax.get_transform('world')
            length = (50*u.au / (415*u.pc)).to(u.deg, u.dimensionless_angles())
            endpoints_x = u.Quantity([0.5*u.arcsec, 0.5*u.arcsec+length]) + leftmost_position
            ax.plot(endpoints_x.to(u.arcsec),
                    ([vrange[0]+2]*2*u.km/u.s).to(u.m/u.s),
                    'r',
                    transform=trans,
                    zorder=100, linewidth=2)
            ax.text(endpoints_x.mean().value,
                    (vrange[0]+3)*scalefactor,
                    "50 au", color='r', transform=trans, ha='center')
            ax.plot(u.Quantity([leftmost_position, rightmost_position]).value,
                    u.Quantity([vcen,vcen]).to(u.m/u.s).value, 'w:', transform=trans)

            origin = offset_to_point(source.ra.deg,
                                     source.dec.deg,
                                     extraction_path)*u.deg

            ax.vlines(origin.to(u.arcsec).value,
                      (vrange[0]-5)*scalefactor,
                      (vrange[1]+5)*scalefactor,
                      color='r', linestyle='--', linewidth=2.0,
                      alpha=0.6, transform=trans)


            ax.set_ylim(ww.wcs_world2pix(0,vrange[0]*scalefactor,0)[1],
                        ww.wcs_world2pix(0,vrange[1]*scalefactor,0)[1])
            ax.set_xlim(good_limits)


            # ax.set_aspect(4)
            ax.set_aspect(2*extracted.data.shape[1]/extracted.data.shape[0])
            #ax.set_aspect('equal')

            ax.coords[1].set_format_unit(u.km/u.s)

            pl.colorbar(im)


            fig.savefig(paths.fpath('pv/{0}/'.format(name, linename) +
                                    basename.replace(".fits",".png")),
                        dpi=200,
                        bbox_inches='tight')

            # overlay a Keplerian velocity curve
            positions = np.linspace(0,1000,200)*u.au
            # this is the 3d velocity, so assumes edge-on
            vel = (((constants.G * 20*u.M_sun)/(positions))**0.5).to(u.m/u.s)
            loc = (positions/(415*u.pc)).to(u.arcsec, u.dimensionless_angles())
            axlims = ax.axis()
            ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=1.0, transform=trans)

            vel = (((constants.G * 10*u.M_sun)/(positions))**0.5).to(u.m/u.s)
            loc = (positions/(415*u.pc)).to(u.arcsec, u.dimensionless_angles())
            vcen = u.Quantity(vcen, u.km/u.s)
            ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=1.0, transform=trans)

            vel = (((constants.G * 5*u.M_sun)/(positions))**0.5).to(u.m/u.s)
            loc = (positions/(415*u.pc)).to(u.arcsec, u.dimensionless_angles())
            vcen = u.Quantity(vcen, u.km/u.s)
            axlims = ax.axis()
            ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'g:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'g:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'g:', linewidth=1.0, alpha=1.0, transform=trans)
            ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'g:', linewidth=1.0, alpha=1.0, transform=trans)

            ax.axis(axlims)

            fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name, linename) +
                                    basename.replace(".fits",".png")),
                        dpi=200,
                        bbox_inches='tight')


            #outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
            #outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
            #extracted = pvextractor.extract_pv_slice(medsub, outflowpath)
            #extracted.writeto('W51e2_PV_outflowaxis_spw{0}.fits'.format(ii), overwrite=True)

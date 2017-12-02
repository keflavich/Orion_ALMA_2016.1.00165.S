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
from lines import lines

# use outflow_meta b/c higher precision than ds9 reg
from line_point_offset import offset_to_point

import pylab as pl

pl.close(1)

diskycoorddict = {}
source = 'sourceI'
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
    ('sourceI', 'sourceI', source, (0,20), 7.0),
   ):

    for linename, linefreq in lines.items():

        diskycoords = diskycoorddict[name]

        fn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout.fits'

        assert 'cutout' in fn
        basename = "{0}_{1}_diskpv.fits".format(source, linename)
        outfn = paths.dpath(os.path.join("pv/", basename))

        cube = (SpectralCube.read(fn)
                .with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                    rest_value=linefreq)
                .mask_out_bad_beams(5))
        cube.allow_huge_operations=True
        cube.beam_threshold = 5000
        med = cube.percentile(25,axis=0)
        medsub = cube - med

        extraction_path = pvextractor.Path(diskycoords, 0.05*u.arcsec)
        extracted = pvextractor.extract_pv_slice(medsub, extraction_path)
        extracted.writeto(outfn, clobber=True)

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
        plotted_slice = (slice(int(plotted_region[1][0]), int(plotted_region[1][1])),
                         slice(None,None),
                        )
        if np.any(np.array(extracted.data[plotted_slice].shape) == 0):
            log.warn("Skipping {0} because it's empty".format(fn))
            continue


        fig = pl.figure(1, figsize=(6,4))
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


        vmin,vmax = (np.nanmin(extracted.data[plotted_slice]),
                     np.nanmax(extracted.data[plotted_slice]))
        im = ax.imshow(extracted.data, cmap='gray_r',
                       vmin=vmin, vmax=vmax*1.1)
        ax.set_xlabel("Offset [\"]")
        ax.set_ylabel("$V_{LSR}$ [km/s]")


        trans = ax.get_transform('world')
        length = (2000*u.au / (5400*u.pc)).to(u.deg, u.dimensionless_angles())
        endpoints_x = u.Quantity([0.5*u.arcsec, 0.5*u.arcsec+length]) + leftmost_position
        ax.plot(endpoints_x.to(u.arcsec),
                ([vrange[0]+2]*2*u.km/u.s).to(u.m/u.s),
                'r',
                transform=trans,
                zorder=100, linewidth=2)
        ax.text(endpoints_x.mean().value,
                (vrange[0]+3)*scalefactor,
                "2000 au", color='r', transform=trans, ha='center')

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
                    bbox_inches='tight')

        # overlay a Keplerian velocity curve
        positions = np.linspace(0,10000,500)*u.au
        # this is the 3d velocity, so assumes edge-on
        vel = (((constants.G * 20*u.M_sun)/(positions))**0.5).to(u.m/u.s)
        loc = (positions/(5410*u.pc)).to(u.arcsec, u.dimensionless_angles())
        vcen = u.Quantity(vcen, u.km/u.s)
        axlims = ax.axis()
        ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'b:', linewidth=1.0, alpha=0.5, transform=trans)

        positions = np.linspace(0,10000,500)*u.au
        vel = (((constants.G * 50*u.M_sun)/(positions))**0.5).to(u.m/u.s)
        loc = (positions/(5410*u.pc)).to(u.arcsec, u.dimensionless_angles())
        vcen = u.Quantity(vcen, u.km/u.s)
        ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=0.5, transform=trans)
        ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(u.m/u.s), 'm:', linewidth=1.0, alpha=0.5, transform=trans)

        ax.axis(axlims)

        fig.savefig(paths.fpath('pv/{0}/keplercurves_'.format(name, linename) +
                                basename.replace(".fits",".png")),
                    bbox_inches='tight')


        #outflow_coords = coordinates.SkyCoord(["19:23:44.127 +14:30:32.30", "19:23:43.822 +14:30:36.64"], unit=(u.hour, u.deg), frame='fk5')
        #outflowpath = pvextractor.Path(outflow_coords, 0.2*u.arcsec)
        #extracted = pvextractor.extract_pv_slice(medsub, outflowpath)
        #extracted.writeto('W51e2_PV_outflowaxis_spw{0}.fits'.format(ii), clobber=True)

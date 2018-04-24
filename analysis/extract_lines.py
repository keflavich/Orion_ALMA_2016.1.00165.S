import numpy as np
import paths
import os
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from spectral_cube import SpectralCube, wcs_utils, tests, Projection
from lines import disk_lines
import radio_beam
from files import b6_hires_cont
from astropy.nddata import Cutout2D
import reproject

import pylab as pl
pl.ioff()

#conthdu = fits.open(paths.dpath('Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy_SourceIcutout.image.tt0.pbcor.fits'))
#conthdu = fits.open(paths.dpath('Orion_SourceI_B6_continuum_r0.5_SourceIcutout.image.tt0.pbcor.fits'))
#conthdu = fits.open(paths.dpath('OrionSourceI_Band6_QA2_continuum_cutout.fits'))
#conthdu = fits.open(paths.dpath('sourceIcutouts/Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4_SourceIcutout.image.tt0.pbcor.fits'))
conthdu = fits.open(paths.dpath(b6_hires_cont))
assert conthdu[0].header['RADESYS'] == 'ICRS'

coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                             unit=(u.hour, u.deg))
mywcs = wcs.WCS(conthdu[0].header).celestial
cutout = Cutout2D(data=conthdu[0].data.squeeze(),
                  wcs=mywcs,
                  position=coord,
                  size=0.5*u.arcsec)

beam = radio_beam.Beam.from_fits_header(conthdu[0].header)
data = cutout.data
new_header = cutout.wcs.to_header()
new_header.update(beam.to_header_keywords())
mywcs = cutout.wcs
conthdu = fits.PrimaryHDU(data=data, header=new_header)



levels = [50, 150, 300, 500]*u.K
levels = [50, 300, 500]*u.K

for robust in (-2, 0.5):
    ftemplate = '/Volumes/external/orion/Orion{1}_only.{2}.robust{robust}.spw{0}.{suffix}_medsub.image.pbcor.fits'


    cont_levels_Jy = levels.to(u.Jy, beam.jtok_equiv(224*u.GHz))
    #print("Levels: {0}".format(([0.001, 0.005, 0.01, 0.02, 0.03, 0.04,
    #                             0.05]*u.Jy).to(u.K, beam.jtok_equiv(224*u.GHz))))

    for band,suffix in (('B3', 'clarkclean10000'),
                        ('B6', 'maskedclarkclean10000')):
        for sourcename in ('SourceI',):# 'BN'):
            for linename, linefreq in disk_lines.items():
                if 'H30a' in linename:
                    vrange = [-120,130]
                elif 'SiO' in linename:
                    vrange = [-50,60]
                else:
                    vrange = [-30,40]

                for spw in (0,1,2,3):
                    filename = ftemplate.format(spw, sourcename, band, suffix=suffix, robust=robust)

                    if os.path.exists(filename):
                        cube = SpectralCube.read(filename)
                        if not cube.wcs.wcs.radesys.lower() == 'icrs':
                            log.exception("File {0} has wrong coordinate system {1}".format(filename, cube.wcs.wcs.radesys))
                            continue
                    else:
                        log.exception("File {0} does not exist".format(filename))
                        continue

                    fmin, fmax = cube.spectral_extrema

                    mywcs = cube.wcs.celestial
                    pixscale = wcs.utils.proj_plane_pixel_area(mywcs)**0.5 * u.deg

                    if linefreq > fmin and linefreq < fmax:


                        mx_fn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.fits').format(linename, sourcename, robust=robust)
                        m0_fn = paths.dpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0.fits').format(linename, sourcename, robust=robust)

                        if not os.path.exists(mx_fn):
                            print("Writing {0} in spw {1}".format(linename, spw))

                            scube = (cube.with_spectral_unit(u.km/u.s,
                                                             velocity_convention='radio',
                                                             rest_value=linefreq)
                                     .spectral_slab(vrange[0]*u.km/u.s,
                                                    vrange[1]*u.km/u.s))

                            scube.beam_threshold = 0.1
                            cubeK = scube.to(u.K)
                            linecubepath = paths.dpath('cubes/Orion{1}_{0}_robust{robust}maskedclarkclean10000_medsub_K.fits').format(linename, sourcename, robust=robust)
                            log.info("Writing {0}".format(linecubepath))
                            cubeK.write(linecubepath, overwrite=True)

                            m0 = cubeK.moment0(axis=0)
                            mx = cubeK.max(axis=0)

                            mx.write(mx_fn,
                                     overwrite=True)
                            m0.write(m0_fn,
                                     overwrite=True)

                            del cubeK
                            del scube
                        else:

                            # the cubes have to exist anyway...
                            mx = Projection.from_hdu(fits.open(mx_fn)[0])
                            m0 = Projection.from_hdu(fits.open(m0_fn)[0])


                        print("Figures for {0} in spw {1}".format(linename, spw))
                        pl.figure(1).clf()

                        mx.quicklook(filename=paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.pdf')
                                     .format(linename, sourcename, robust=robust), aplpy_kwargs={'figure':pl.figure(1)})
                        mx.FITSFigure.show_grayscale(invert=True)
                        mx.FITSFigure.show_contour(conthdu, levels=cont_levels_Jy,
                                                   colors=['r']*10, linewidths=0.75)
                        mx.FITSFigure.colorbar.set_axis_label_text("$T_B$ [K]")
                        mx.FITSFigure.save(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak.pdf')
                                           .format(linename, sourcename, robust=robust))

                        mx.FITSFigure.add_beam()
                        mx.FITSFigure.beam.set_facecolor('none')
                        mx.FITSFigure.beam.set_edgecolor('b')

                        pl.figure(1).clf()

                        contdata,_ = reproject.reproject_interp(conthdu, mx.hdu.header)

                        extent = [-mx.shape[1]/2*pixscale.to(u.arcsec).value,
                                  mx.shape[1]/2*pixscale.to(u.arcsec).value,
                                  -mx.shape[0]/2*pixscale.to(u.arcsec).value,
                                  mx.shape[0]/2*pixscale.to(u.arcsec).value,]

                        ax = pl.figure(1).gca()
                        im = ax.imshow(mx.value, cmap='gray_r',
                                       interpolation='none', origin='lower',
                                       extent=extent,
                                       zorder=0,
                                      )

                        ebm = beam.ellipse_to_plot(extent[0]+0.05, extent[2]+0.05, 1./3600*u.deg)
                        ebm.set_facecolor('none')
                        ebm.set_edgecolor('b')
                        ax.add_patch(ebm)

                        con = ax.contour(contdata, levels=cont_levels_Jy.value,
                                         colors=['r']*10, extent=extent, linewidths=0.75,
                                         zorder=50,
                                        )
                        cb = pl.colorbar(mappable=im)
                        cb.set_label("$T_B$ [K]")
                        ax.set_xlabel("RA offset (\")")
                        ax.set_ylabel("Dec offset (\")")

                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak_offset.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))
                        if mx.value.max() < 250:
                            line_levels = [50,100,150,200]
                        elif mx.value.max() > 1000:
                            line_levels = [500, 750, 1000, 1250]
                        else:
                            line_levels = [200, 400, 600, 800]
                        con2 = ax.contour(mx.value, levels=line_levels,
                                          zorder=25,
                                          colors=['w']*5, extent=extent, linewidths=0.75)
                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak_offset_contours.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))
                        con.set_alpha(0)
                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_peak_offset_contours_nocont.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))

                        pl.figure(1).clf()

                        m0.quicklook(filename=paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0.pdf')
                                     .format(linename, sourcename, robust=robust), aplpy_kwargs={'figure':pl.figure(1)})
                        m0.FITSFigure.show_grayscale(invert=True)
                        m0.FITSFigure.add_beam()
                        m0.FITSFigure.beam.set_facecolor('none')
                        m0.FITSFigure.beam.set_edgecolor('b')

                        m0.FITSFigure.show_contour(conthdu, levels=cont_levels_Jy,
                                                   colors=['r']*10, linewidths=0.75)
                        m0.FITSFigure.colorbar.set_axis_label_text("$\int T_B \mathrm{d}v$ [K km s$^{-1}$]")
                        m0.FITSFigure.save(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0.pdf')
                                           .format(linename, sourcename, robust=robust))

                        pl.figure(1).clf()


                        ax = pl.figure(1).gca()
                        im = ax.imshow(m0.value, cmap='gray_r',
                                       interpolation='none', origin='lower',
                                       extent=extent,
                                       zorder=0,
                                      )

                        ebm = beam.ellipse_to_plot(extent[0]+0.05, extent[2]+0.05, 1./3600*u.deg)
                        ebm.set_facecolor('none')
                        ebm.set_edgecolor('b')
                        ax.add_patch(ebm)

                        con = ax.contour(contdata, levels=cont_levels_Jy.value,
                                         zorder=50,
                                         colors=['r']*10, extent=extent, linewidths=0.75)
                        cb = pl.colorbar(mappable=im)
                        cb.set_label("$\int T_B \mathrm{d}v$ [K km s$^{-1}$]")
                        ax.set_xlabel("RA offset (\")")
                        ax.set_ylabel("Dec offset (\")")

                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0_offset.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))

                        plotted_center = mywcs.wcs_pix2world(m0.shape[1]/2, m0.shape[0]/2, 0)
                        plotted_center_coord = coordinates.SkyCoord(*plotted_center, unit=(u.deg,
                                                                                           u.deg),
                                                                    frame=mywcs.wcs.radesys.lower())
                        print("Center position of {2} image is: {0}  {1}"
                              .format(plotted_center_coord.to_string('hmsdms'),
                                      plotted_center_coord.to_string('hmsdms', sep=':'),
                                      band
                                     )
                             )

                        stderr = stats.mad_std(m0.value, ignore_nan=True)
                        con.set_alpha(0)
                        con2 = ax.contour(m0.value, levels=np.array([5,10,15,20,25])*stderr,
                                          zorder=25,
                                          colors=['w']*5, extent=extent, linewidths=0.75)
                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0_offset_contours_nocont.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))

                        con2.set_alpha(1)
                        con.set_alpha(1)
                        for coll in con2.collections:
                            coll.set_color('w')
                        pl.savefig(paths.fpath('moments/Orion{1}_{0}_robust{robust}.maskedclarkclean10000_medsub_K_moment0_offset_contours.pdf')
                                   .format(linename, sourcename,
                                           robust=robust))




                        pl.close('all')


                        del mx, m0
                    del cube


#    beam = radio_beam.Beam.from_fits_header(conthdu[0].header)
#    print("Levels: {0}".format(([0.001, 0.005, 0.01, 0.02, 0.03, 0.04,
#                                 0.05]*u.Jy).to(u.K, beam.jtok_equiv(224*u.GHz))))

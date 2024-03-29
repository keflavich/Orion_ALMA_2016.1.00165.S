import numpy as np
import regions
from astropy import units as u
from astropy.modeling import models, fitting
from astropy import wcs
from astropy import coordinates
from astropy.io import fits
from astropy import convolution
from spectral_cube import SpectralCube
import pvextractor
from gaussfit_catalog.core import gaussfit_image
import paths
import os
import pylab as pl
import shapely.geometry as geom
from show_pv import show_keplercurves, show_pv
import imp
import edge_on_ring_velocity_model
from astropy.utils.console import ProgressBar
from constants import vcen as assumed_vcen, origin, d_orion
import tracemalloc
from astropy import log

imp.reload(edge_on_ring_velocity_model)
from edge_on_ring_velocity_model import thindiskcurve, thindiskcurve_fitter, trace_plotter

first_snapshot = tracemalloc.take_snapshot()

if 'cached_gaussfit_results' not in locals():
    cached_gaussfit_results = {}

fitresult_list = []

for linename,(vmin,vmax),limits,(cenx, ceny) in (
    #('Unknown_4', (-15, 27), (-0.1, 0.1, -0.12, 0.12), (65.1, 60.5)),
    #('H2Ov2=1_5(5,0)-6(4,3)', (-28, 38), (-0.2, 0.2, -0.2, 0.2), (68.8, 65.5)),
    #('Unknown_1', (-15, 27), (-0.1, 0.1, -0.12, 0.12), (65.1, 60.5)),
    #('Unknown_2', (-15, 27), (-0.1, 0.1, -0.12, 0.12), (65.1, 60.5)),
    ('SiOv=1_5-4', (-30, 45), (-0.2, 0.2, -0.2, 0.2), (65.1, 60.5)),
   ):

    regs = regions.Regions.read(paths.rpath('velo_centroid_guesses_{linename}.reg').format(linename=linename))

    guesses = {}

    for reg in regs:
        vel = float(reg.meta['text'].strip("{}")
                    if 'text' in reg.meta
                    else reg.meta['label']
                   )
        if vel in guesses:
            guesses[vel].append(reg)
        else:
            guesses[vel] = [reg]

    velocities = np.array(sorted(guesses.keys()))*u.km/u.s

    cubefn = paths.dpath('cubes/OrionSourceI_{linename}_robust0.5maskedclarkclean10000_medsub_K.fits'
                         .format(linename=linename))
    basename = os.path.splitext(os.path.basename(cubefn))[0]
    cube = SpectralCube.read(cubefn)

    # hard-coded icrs radecsys assumed
    assert cube.header['RADESYS'] == 'ICRS'

    vdiff = np.abs(np.diff(cube.spectral_axis).mean())
    rms = cube.std()
    weights = np.ones(cube.shape[1:]) * rms.value**2

    try:
        beam = cube.beam
    except AttributeError:
        beam = cube.average_beams(1)
    pixscale = wcs.utils.proj_plane_pixel_area(cube.wcs.celestial)**0.5 * u.deg
    bmmaj_px = (beam.major / pixscale).decompose().value
    bmmin_px = (beam.minor / pixscale).decompose().value
    max_radius_in_beams = 1.25
    max_offset_in_beams = 1.5

    STDDEV_TO_FWHM = np.sqrt(8*np.log(2))

    results = {}

    if linename in cached_gaussfit_results:
        print("Loading {0} from cache".format(linename))
        results = cached_gaussfit_results[linename]
    else:
        print("Fitting {0}, which is not in cache".format(linename))
        pb = ProgressBar(cube.shape[0])

        for vel,vslice in (zip(cube.spectral_axis, cube)):
            closest = np.argmin(np.abs(vel-velocities))
            if np.abs(velocities[closest] - vel) > vdiff:
                #print("Skipping velocity {0}, closest is {1} -> {2}".format(vel, closest,
                #                                                            velocities[closest]))
                pb.update()
                continue

            thisvel = velocities[closest].value
            guess_regs = guesses[thisvel]

            ampguess = vslice.max().value

            model_list = []
            for reg in guess_regs:

                p_init = models.Gaussian2D(amplitude=ampguess,
                                           x_mean=reg.center.x,
                                           y_mean=reg.center.y,
                                           x_stddev=bmmaj_px/STDDEV_TO_FWHM*0.75,
                                           y_stddev=bmmin_px/STDDEV_TO_FWHM*0.75,
                                           theta=beam.pa.to(u.deg).value,
                                           bounds={'x_stddev':(bmmin_px/STDDEV_TO_FWHM*0.5,
                                                               bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                                   'y_stddev':(bmmin_px/STDDEV_TO_FWHM*0.5,
                                                               bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                                   'x_mean':(reg.center.x-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                             reg.center.x+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                                   'y_mean':(reg.center.y-max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM,
                                                             reg.center.y+max_offset_in_beams*bmmaj_px/STDDEV_TO_FWHM),
                                                   'amplitude':(ampguess*0.1, ampguess*2.1)
                                                  }
                                          )
                assert not p_init._has_units
                model_list.append(p_init)

            composite_model = model_list[0]
            if len(model_list) > 1:
                for mod in model_list[1:]:
                    composite_model += mod

            fit_result = gaussfit_image(vslice.value, composite_model,
                                        weights=weights, plot=True)
            results[thisvel] = fit_result

            if not os.path.exists(paths.fpath('velcentroid/diagnostics/{0}/'.format(basename))):
                os.mkdir(paths.fpath('velcentroid/diagnostics/{0}/').format(basename))
            pl.savefig(paths.fpath('velcentroid/diagnostics/{0}/{1}.png')
                       .format(basename, thisvel))

            pb.update()

        cached_gaussfit_results[linename] = results


    diskend_regs = regions.Regions.read(paths.rpath('diskends.reg'))
    diskends = coordinates.SkyCoord([reg.center for reg in diskend_regs]).icrs

    #center = coordinates.SkyCoord(diskends.ra.mean(), diskends.dec.mean(),
    #                              frame=diskends.frame)
    center_cont = regions.Regions.read(paths.rpath('sourceI_center.reg'))[0].center.transform_to(coordinates.ICRS)
    center_U = regions.Regions.read(paths.rpath('sourceI_center.reg'))[1].center.transform_to(coordinates.ICRS)
    center_h2o = regions.Regions.read(paths.rpath('sourceI_center.reg'))[2].center.transform_to(coordinates.ICRS)

    center = center_cont

    ref_cen_x, ref_cen_y = cube.wcs.celestial.wcs_world2pix(center.ra.deg,
                                                            center.dec.deg, 0)
    # had to switch to pixel centering because coordinate transformations
    # aren't accurate enough (and/or, don't match ds9)
    ref_cen_x = cenx
    ref_cen_y = ceny

    #assert np.abs(ref_cen_x - 66) < 2
    #assert np.abs(ref_cen_y - 63) < 2

    def offset_to_point(xx, yy):
        line = geom.LineString(zip(diskends.ra.deg, diskends.dec.deg))
        point = geom.Point(xx, yy)
        return line.project(point)

    ref_offset = offset_to_point(center.ra.deg, center.dec.deg)

    # Add the data as a background image
    pvfn = 'pv/sourceI_{0}_B6_robust0.5_diskpv_0.01.fits'.format(linename)
    fh = fits.open(paths.dpath(pvfn))
    data = fh[0].data
    header = fh[0].header

    ww = wcs.WCS(fh[0].header)

    ww.wcs.cdelt[1] /= 1000.0
    ww.wcs.crval[1] /= 1000.0
    ww.wcs.cunit[1] = ' '#u.km/u.s
    ww.wcs.cdelt[0] *= 3600
    ww.wcs.cunit[0] = u.arcsec
    ww.wcs.crval[0] = -origin.to(u.arcsec).value

    # #!@#$!@#$@!%@#${^(@#$)%#$(
    ww.wcs.set()


    snap3 = tracemalloc.take_snapshot()
    diff = snap3.compare_to(first_snapshot, 'lineno')
    mem_used = sum([dd.size_diff for dd in diff])*u.B
    log.info(f"Used {mem_used.to(u.GB)} for line {linename} prior to plotting")

    fig1 = pl.figure(1)
    pl.clf()
    ax1 = fig1.add_subplot(1,2,1)
    pvwcs = pvextractor.utils.wcs_slicing.slice_wcs(cube.wcs, pixscale)
    pvwcs.wcs.cdelt[0] *= 3600
    pvwcs.wcs.crval[0] = -origin.to(u.arcsec).value
    pvwcs.wcs.cunit[0] = u.arcsec
    pvwcs.wcs.crval[1] /= 1e3
    # not clear why this is needed?
    pvwcs.wcs.crval[1] += assumed_vcen.to(u.km/u.s).value
    pvwcs.wcs.cdelt[1] /= 1e3
    pvwcs.wcs.cunit[1] = ' ' # force wcs.set to NOT revert to m/s
    print('pvwcs: cdelt={0} crval={1} crpix={2}'.format(pvwcs.wcs.cdelt,
                                                        pvwcs.wcs.crval,
                                                        pvwcs.wcs.crpix))
    ax2 = fig1.add_subplot(1,2,2, projection=ww)

    trans = ax2.get_transform('world')

    diskends_x, diskends_y = cube.wcs.celestial.wcs_world2pix(diskends.ra.deg,
                                                              diskends.dec.deg,
                                                              0)
    ax1.plot((diskends_x-ref_cen_x)*pixscale.to(u.arcsec),
             (diskends_y-ref_cen_y)*pixscale.to(u.arcsec), 'k-',
             linewidth=3, alpha=0.2, zorder=-100)
    ax1.plot(0, 0, 'ko', markersize=5, alpha=0.5, zorder=-50)

    loc = pl.matplotlib.ticker.MultipleLocator(base=0.04)
    ax1.xaxis.set_major_locator(loc)

    cmap = pl.cm.Spectral_r
    cmap = pl.cm.Spectral
    norm = pl.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    offset_fits_arcsec = []
    vels = np.sort(np.array([v for v in results]))

    gfit_points = []

    for vel in vels:
        fit_result = results[vel]
        fitted = fit_result[0]
        n_submodels = fitted.n_submodels()
        color = cmap(norm(vel))
        if n_submodels > 1:
            offsets_ = []
            for ii in range(n_submodels):
                xcen, ycen = (getattr(fitted, 'x_mean_{0}'.format(ii)),
                              getattr(fitted, 'y_mean_{0}'.format(ii)))
                ax1.plot((xcen-ref_cen_x)*pixscale.to(u.arcsec),
                         (ycen-ref_cen_y)*pixscale.to(u.arcsec),
                         color=color,
                         marker='o')
                ra,dec = cube.wcs.celestial.wcs_pix2world(xcen, ycen, 0)
                offset = offset_to_point(ra, dec)
                gfit_points.extend(ax2.plot((offset-ref_offset)*3600, vel,
                                            color=color, marker='s',
                                            markersize=5,
                                            markeredgecolor='k',
                                            transform=trans))
                offsets_.append((offset-ref_offset)*3600)
            offset_fits_arcsec.append(np.mean(offsets_))
        else:
            xcen, ycen = fitted.x_mean, fitted.y_mean
            ax1.plot((xcen-ref_cen_x)*pixscale.to(u.arcsec),
                     (ycen-ref_cen_y)*pixscale.to(u.arcsec), color=color, marker='o')
            ra,dec = cube.wcs.celestial.wcs_pix2world(xcen, ycen, 0)
            offset = offset_to_point(ra, dec)
            gfit_points.extend(ax2.plot((offset-ref_offset)*3600, vel,
                                        color=color, marker='s',
                                        markersize=5,
                                        markeredgecolor='k',
                                        transform=trans)
                              )
            offset_fits_arcsec.append((offset-ref_offset)*3600)


    show_keplercurves(ax2, 0*u.deg, 150, assumed_vcen, yaxis_unit=u.km/u.s, radii={},
                      masses=[15, ], linestyles=[':',], colors=['r'],
                      trans=ax2.get_transform('world'),
                     )

    # xx_thindisk, yy_thindisk = thindiskcurve(mass=20*u.M_sun, rmin=30*u.au, rmax=80*u.au)
    # ax2.plot((xx_thindisk / d_orion).to(u.arcsec, u.dimensionless_angles()),
    #          yy_thindisk + assumed_vcen,
    #          'k:',
    #          transform=trans)
    xx_thindisk, yy_thindisk = thindiskcurve(mass=15.5*u.M_sun, rmin=17*u.au,
                                             rmax=66*u.au)
    thindiskline = ax2.plot((xx_thindisk / d_orion).to(u.arcsec,
                                                       u.dimensionless_angles()),
                            yy_thindisk + assumed_vcen, 'k-', transform=trans)

    xlim, ylim = pvwcs.wcs_world2pix([-0.2,0.2], [vmin-3, vmax+3], 0)
    xlim, ylim = ww.wcs_world2pix([-0.2,0.2], [vmin-3, vmax+3], 0)
    ax2.set_ylim(*ylim)
    ax2.set_xlim(*xlim)
    ax2.yaxis.tick_right() # mpl version; incompatible with wcsaxes
    ax2.coords[1].set_ticklabel_position('r')
    ax2.coords[1].set_axislabel_position('r')
    ax2.set_ylabel("$V_{LSR}$ [km s$^{-1}$]")
    ax2.set_xlabel("Offset Position (arcsec)")
    ax1.set_xlabel("Offset RA (arcsec)")
    ax1.set_ylabel("Offset Dec (arcsec)")
    ax1.axis(limits)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(90)
    for tick in ax2.get_xticklabels():
        tick.set_rotation(90)

    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots.pdf'.format(linename=linename)),
               bbox_inches='tight')

    im = ax2.imshow(data, cmap='gray_r', interpolation='none', vmax=0.02,
                    vmin=-0.0001, origin='lower',
                    #transform=ax2.get_transform('world'),
                   )
    ax2.set_aspect('auto')
    ax2.set_ylim(*ylim)
    ax2.set_xlim(*xlim)

    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_withbackground.pdf'.format(linename=linename)),
               bbox_inches='tight')

    #fig3,ax3 = show_pv(data, ww,
    #                   origin, vrange=[-40,55], vcen=assumed_vcen,
    #                   imvmin=-0.001, imvmax=0.020)

    #trans3 = ax3.get_transform('world')
    #overlay_thindiskline = ax3.plot((xx_thindisk / d_orion).to(u.arcsec,
    #                                                           u.dimensionless_angles()),
    #                                 (yy_thindisk + assumed_vcen).to(u.m/u.s), 'k-', transform=trans3)


    # previously had an
    #if 'Unknown' in linename:
    # statement here because these are the only ones that represent disks,
    # but since we're using the *average* offset position, maybe this is OK?

    offsets_au = (u.Quantity(offset_fits_arcsec, u.arcsec) *
                  d_orion).to(u.au, u.dimensionless_angles())
    print("Beginning thin disk fitting for {0}".format(linename))
    if False:
        fitresult, ci, trace, minimizer = thindiskcurve_fitter(
            xsep=np.array(offsets_au), velo=vels, mguess=15.5*u.M_sun, rinner=17,
            router=66, conf_int=True, npix=500)
        fitresult_list.append(fitresult)
        trace_plotter(ci, trace, fitresult, minimizer)
    else:
        fitresult = thindiskcurve_fitter(xsep=np.array(offsets_au), velo=vels,
                                         mguess=15.5*u.M_sun, rinner=17,
                                         router=66, conf_int=False, npix=500)
        fitresult_list.append(fitresult)

    for line in thindiskline:
        line.set_visible(False)

    xx_thindisk, yy_thindisk = thindiskcurve(mass=fitresult.params['mass']*u.M_sun,
                                             rmin=fitresult.params['rinner']*u.au,
                                             rmax=fitresult.params['router']*u.au,
                                            )


    lines = ax2.plot((xx_thindisk / d_orion).to(u.arcsec, u.dimensionless_angles()),
                     yy_thindisk + assumed_vcen,
                     'k-',
                     transform=trans,
                     label=('$M={0:0.1f}$\n$R_{{in}}={1:d}$\n$R_{{out}}={2:d}$'
                            .format(fitresult.params['mass'].value,
                                    int(fitresult.params['rinner'].value),
                                    int(fitresult.params['router'].value),
                                   )
                           )
                    )

    xmin,xmax = ax2.get_xlim()
    ax2.hlines(fitresult.params['vcen'].value, xmin, xmax, linestyle='--', color='k',
               alpha=0.5, zorder=500, transform=trans)

    leg = pl.legend(loc='upper left', fontsize=12, handlelength=1.0)

    im.set_visible(False)

    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel.pdf'.format(linename=linename)),
               bbox_inches='tight')

    # show the average fitted positions
    avgposns = ax2.plot(offset_fits_arcsec, vels, marker='o', markersize=5,
                        linestyle='none', markerfacecolor='k', alpha=0.5,
                        markeredgecolor='w', transform=trans)

    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_withavgs.pdf'.format(linename=linename)),
               bbox_inches='tight')

    im.set_visible(True)
    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_withavgs_withbackground.pdf'.format(linename=linename)),
               bbox_inches='tight')
    im.set_visible(False)


    tracemalloc.start()

    snap1 = tracemalloc.take_snapshot()

    for mass in (15.5, 19, 5, 10, 15, 20, ):
        for line in lines:
            line.remove()
            #line.set_visible(False)
        leg.remove()
        del leg, lines

        fitresult = thindiskcurve_fitter(xsep=np.array(offsets_au),
                                         velo=vels,
                                         mguess=mass*u.M_sun,
                                         rinner=17,
                                         router=66,
                                         fixedmass=True,
                                        )
        assert fitresult.params['mass'].value == mass

        snap2 = tracemalloc.take_snapshot()
        diff = snap2.compare_to(snap1, 'lineno')
        mem_used = sum([dd.size_diff for dd in diff])*u.B
        log.info(f"Used {mem_used.to(u.GB)} with mass={mass} after fitresult")


        xx_thindisk, yy_thindisk = thindiskcurve(mass=fitresult.params['mass']*u.M_sun,
                                                 rmin=fitresult.params['rinner']*u.au,
                                                 rmax=fitresult.params['router']*u.au,
                                                )

        lines = ax2.plot((xx_thindisk / d_orion).to(u.arcsec, u.dimensionless_angles()),
                         yy_thindisk + assumed_vcen,
                         'k-',
                         transform=trans,
                         label=('$M={0:0.1f}$\n$R_{{in}}={1:d}$\n$R_{{out}}={2:d}$'
                                .format(fitresult.params['mass'].value,
                                        int(fitresult.params['rinner'].value),
                                        int(fitresult.params['router'].value),
                                       )
                               )
                        )
        leg = pl.legend(lines,
                        [line.get_label() for line in lines],
                        loc='upper left',
                        fontsize=12, handlelength=1.0)
        pl.tight_layout()
        pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_{mass}msun_withavgs.pdf'
                               .format(mass=mass, linename=linename)),
                   bbox_inches='tight'
                  )


        del fitresult, xx_thindisk, yy_thindisk

        snap2 = tracemalloc.take_snapshot()
        diff = snap2.compare_to(snap1, 'lineno')
        mem_used = sum([dd.size_diff for dd in diff])*u.B
        log.info(f"Used {mem_used.to(u.GB)} with mass={mass}")

    im.set_visible(True)
    for line in lines:
        line.set_visible(False)
    leg.remove()

    xx,vv,pvd = thindiskcurve(mass=15*u.M_sun,
                              rmin=25*u.au,
                              rmax=65*u.au,
                              pvd=True,
                             )
    xx_as = (xx/d_orion).to(u.arcsec, u.dimensionless_angles())

    conv_beam = convolution.Gaussian2DKernel(0.04*u.arcsec/np.diff(xx_as).mean(),
                                             1.5*u.km/u.s/np.diff(vv).mean(),
                                            )
    sm_pvd = convolution.convolve_fft(pvd, conv_beam)


    ax1.cla()
    #ax1.remove()
    #ax1 = fig1.add_subplot(1,2,1, projection=ww)
    ax1.imshow(sm_pvd, cmap='gray_r', extent=[xx_as.value.min(),
                                              xx_as.value.max(),
                                              (vv+assumed_vcen).value.min(),
                                              (vv+assumed_vcen).value.max()],
               interpolation='none', origin='lower',
               vmin=0,
               vmax=sm_pvd[np.abs(vv)<20*u.km/u.s].max(),
               label=('$M={0:0.1f}$\n$R_{{in}}={1:d}$\n$R_{{out}}={2:d}$'
                      .format(15, 25, 65))
              )

    keplines = show_keplercurves(ax1, 0*u.deg, 150, assumed_vcen,
                                 yaxis_unit=u.km/u.s, radii={}, masses=[15, ],
                                 linestyles=[':',], colors=['r'],
                                 trans=ax1.transData,
                     )

    xlim_as, ylim_kms = [-0.2,0.2], [vmin-3, vmax+3]
    ax1.set_xlim(*xlim_as)
    ax1.set_ylim(*ylim_kms)

    ax1.set_aspect('auto')
    #ax1.legend(('$M={0:0.1f}$\n$R_{{in}}={1:d}$\n$R_{{out}}={2:d}$'
    #                  .format(15, 25, 65)),
    #           [im],
    #           loc='best')

    loc = pl.matplotlib.ticker.MultipleLocator(base=0.1)
    ax1.xaxis.set_major_locator(loc)

    ax1.set_ylabel("$V_{LSR}$ [km s$^{-1}$]")
    ax1.set_xlabel("Offset Position (arcsec)")


    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_{mass}msun_withavgs_comparepv.pdf'
                           .format(mass=15, linename=linename)),
               bbox_inches='tight')

    for feature in gfit_points + avgposns:
        feature.set_visible(False)

    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_{mass}msun_comparepv.pdf'
                           .format(mass=15, linename=linename)),
               bbox_inches='tight')
    pl.savefig(paths.fpath('velcentroid/{linename}_pp_pv_plots_fittedmodel_{mass}msun_comparepv.png'
                           .format(mass=15, linename=linename)),
               bbox_inches='tight')


    fig1.clf()

    snap3 = tracemalloc.take_snapshot()
    diff = snap3.compare_to(first_snapshot, 'lineno')
    mem_used = sum([dd.size_diff for dd in diff])*u.B
    log.info(f"Used {mem_used.to(u.GB)} for line {linename}")

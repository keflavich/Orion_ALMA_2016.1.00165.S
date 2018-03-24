import numpy as np
import paths
import pvextractor
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
import regions
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy import modeling
from astropy import constants
from scipy import optimize
import imp
import edge_on_ring_velocity_model
from constants import d_orion, vcen, origin
import show_pv
import pylab as pl

imp.reload(edge_on_ring_velocity_model)
imp.reload(show_pv)


#source = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
#                             unit=(u.hour, u.deg))
# fitted *disk* center from diskmodel
all_voff_results = {}

for fn, vmin, vmax, savename, rms, radii, start_dv in [
                                             #('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_robust0.5_diskpv.fits', -0.0005, 0.055,
                                             # 'H2O_kepler_SeifriedPlot.png', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust0.5_diskpv_0.01.fits', -0.0005, 0.025,
                                              'H2O_kepler_SeifriedPlot_0.01arcsec.pdf', 1*u.mJy, [10,100], 7.5),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust0.5_diskpv_0.1.fits', -0.0005, 0.055,
                                              'H2O_kepler_SeifriedPlot_0.1arcsec.pdf', 1*u.mJy, [10,100], 5.0),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.1.fits', -0.0005, 0.03,
                                              'H2O_kepler_SeifriedPlot_0.1arcsec_robust-2.pdf', 0.7*u.mJy, [10,100], 5.0),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.2.fits', -0.0005, 0.03,
                                              'H2O_kepler_SeifriedPlot_0.2arcsec_robust-2.pdf', 0.7*u.mJy, [10,100], 5.0),
                                             #('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.01.fits', -0.0005, 0.048,
                                             # 'H2O_kepler_SeifriedPlot_0.01arcsec.pdf', 1*u.mJy, [10,100], 5.0),
                                             ('pv/sourceI_29SiOv=0_5-4_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              '29SiOv0_5-4_kepler_SeifriedPlot.pdf', 2*u.mJy, [10,100], 5.0),
                                             ('pv/sourceI_SiS_12-11_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              'SiS_12-11_kepler_SeifriedPlot.pdf', 1*u.mJy, [30,200], 5.0),
                                             ('pv/sourceI_Unknown_1_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_1_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80], 5.0),
                                             ('pv/sourceI_Unknown_4_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_4_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80], 5.0),
                                             ('pv/sourceI_Unknown_5_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_5_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80], 5.0),
                                             ('pv/sourceI_29SiOv=0_2-1_B3_robust-2_diskpv_0.01.fits', -0.05, 1,
                                              '29SiOv0_2-1_kepler_SeifriedPlot.pdf', 1*u.mJy, [10,100], 5.0),
                                             ('pv/sourceI_SiOv=1_5-4_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              'SiOv1_5-4_kepler_SeifriedPlot.pdf', 2*u.mJy, [10,100], 5.0),
                                            ]:
    print(fn, vmin, vmax, savename, rms)
    fh = fits.open(paths.dpath(fn))
    data = fh[0].data
    header = fh[0].header

    ww = wcs.WCS(fh[0].header)

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

    vcen_ypix = int(np.round(ww.sub([2]).wcs_world2pix([vcen.value], 0)[0][0]))
    vcen_ypix_left = int(np.round(ww.sub([2]).wcs_world2pix([(vcen-start_dv*u.km/u.s).to(u.m/u.s).value], 0)[0][0]))
    vcen_ypix_right = int(np.round(ww.sub([2]).wcs_world2pix([(vcen+start_dv*u.km/u.s).to(u.m/u.s).value], 0)[0][0]))
    xcen_xpix = int(np.round(ww.sub([1]).wcs_world2pix([0], 0)[0][0]))

    # rms ~ 1 mJy
    threshold = rms * 5

    def descend(xpix, ypix=vcen_ypix, direction=1, threshold=threshold):
        value = u.Quantity(data[ypix, xpix], u.Jy)
        if not np.isfinite(value):
            return np.nan
        if value < threshold:
            #print(xpix, ypix, direction, value, threshold)
            return np.nan
        while value > threshold:
            ypix += direction
            value = u.Quantity(data[ypix, xpix], u.Jy)
            if ypix < 0:
                break

        return ypix

    cdelt_sign = int(np.sign(ww.wcs.cdelt[1]))

    xoffs_as = u.Quantity(ww.sub([1]).wcs_pix2world(np.arange(data.shape[1]), 0), u.arcsec).squeeze()
    xoffs_au = (xoffs_as*415*u.pc).to(u.au, u.dimensionless_angles())
    all_voff_results[savename] = {'xoffs_as': xoffs_as, 'xoffs_au': xoffs_au}

    for threshold_ in (3, 7, 5):
        threshold = threshold_ * rms
        voffs = u.Quantity(ww.sub([2]).wcs_pix2world([descend(ii,
                                                              direction=-cdelt_sign if ii < xcen_xpix else cdelt_sign,
                                                              ypix=vcen_ypix_left if ii < xcen_xpix else vcen_ypix_right,
                                                              threshold=threshold,
                                                             )
                                                      for ii in range(data.shape[1])], 0),
                           u.m/u.s).squeeze()
        # store individual threshold for later comparison
        all_voff_results[savename]['voffs_{0}'.format(threshold_)] = voffs

    all_voff_results[savename]['voffs'] = voffs

    fig,ax = show_pv.show_pv(data, ww,
                             origin, vrange=[-40,55], vcen=vcen,
                             imvmin=vmin, imvmax=vmax)

    ax.plot(xoffs_as, voffs, 'o-', transform=ax.get_transform('world'), markersize=3, markeredgecolor='b',
            zorder=200, alpha=0.9)
    maxdist=150*u.au
    lines = show_pv.show_keplercurves(ax, origin, maxdist, vcen, masses=[15, ],
                                      linestyles=['-'], colors=['r'],
                                      radii={15: (radii, ('m','m'))})
    ax.set_aspect(2)

    x1,x2 = ww.sub([1]).wcs_world2pix([-0.30,0.30],0)[0]
    y1,y2 = ww.sub([2]).wcs_world2pix([-45*scalefactor,55*scalefactor],0)[0]
    ax.axis([x1,x2,y1,y2])


    fig.savefig(paths.fpath('pv/{0}'.format(savename)),
                dpi=200,
                bbox_inches='tight')

    ax.plot(xoffs_as, all_voff_results[savename]['voffs_3'], 's-', transform=ax.get_transform('world'), markersize=3, markeredgecolor='r',
            markerfacecolor='r',
            zorder=250, alpha=0.9)
    ax.plot(xoffs_as, all_voff_results[savename]['voffs_7'], 'h-', transform=ax.get_transform('world'), markersize=3, markeredgecolor='g',
            markerfacecolor='g',
            zorder=250, alpha=0.9)
    fig.savefig(paths.fpath('pv/{0}'.format(savename.replace(".pdf","_threshold_comparison.pdf"))),
                dpi=200,
                bbox_inches='tight')

    #lines = show_pv.show_keplercurves(ax, origin, maxdist, vcen, masses=[15, ],
    #                                  linestyles=['-'], colors=['g'],
    #                                  radii={},
    #                                  show_other_powerlaws=True
    #                                 )
    #fig.savefig(paths.fpath('pv/{0}'.format(savename.replace(".pdf","_withalpha1.pdf"))),
    #            dpi=200,
    #            bbox_inches='tight')


    for line in ax.get_lines() + ax.collections:
        line.set_visible(False)

    xpv, ypv, pv_15msun = edge_on_ring_velocity_model.thindiskcurve(mass=15*u.M_sun,
                                                                    rmin=20*u.au,
                                                                    rmax=70*u.au,
                                                                    vgrid=np.linspace(-35,
                                                                                      35,
                                                                                      200)*u.km/u.s,
                                                                    pvd=True)

    #pv_15msun_c = convolve_fft(pv_15msun, Gaussian2DKernel(0.5))
    xpv_as = (xpv / d_orion).to(u.arcsec, u.dimensionless_angles())
    ax.contour(xpv_as, (ypv+vcen).to(u.m/u.s), pv_15msun, levels=[1, 25],
               colors=['r','b'], transform=ax.get_transform('world')
              )


    EODsavename = savename.replace('Seifried', 'EdgeOnDiskPV')
    fig.savefig(paths.fpath('pv/{0}'.format(EODsavename)),
                dpi=200,
                bbox_inches='tight')


# Post-referee statistical fits to the envelopes

for key,xmin in [('SiS_12-11_kepler_SeifriedPlot.pdf', 25*u.au),
                 ('H2O_kepler_SeifriedPlot_0.01arcsec.pdf', 13*u.au)]:
    fitter = modeling.fitting.LevMarLSQFitter()

    xoffs_au = all_voff_results[key]['xoffs_au']
    voffs = all_voff_results[key]['voffs']

    def vcen_resid(vcen):
        vcen = u.Quantity(vcen, u.km/u.s)
        # first, interpolate across the NaNs to avoid bad values
        interped_voffs = np.interp(xoffs_au, xoffs_au[np.isfinite(voffs)], voffs[np.isfinite(voffs)])
        # then, interpolate the negative offsets onto the positive grid
        voffs_neg = np.interp(xoffs_au[xoffs_au>xmin],
                              -xoffs_au[xoffs_au<-xmin][::-1],
                              interped_voffs[xoffs_au<-xmin][::-1]) * voffs.unit
        voffs_pos = voffs[xoffs_au>xmin]

        # compute the residual
        resid = (voffs_pos[np.isfinite(voffs_pos)] - vcen) - (vcen - voffs_neg[np.isfinite(voffs_pos)])

        return np.sum(resid**2).value

    optimal_vcen = u.Quantity(optimize.fmin(vcen_resid, 5.5*u.km/u.s)[0], u.km/u.s)
    print("optimal vcen for {key} is {optimal_vcen}".format(**locals()))


    xoffs_au_tofit = np.abs(xoffs_au)
    voffs_tofit = np.abs(voffs-optimal_vcen)

    ok = (xoffs_au_tofit > xmin) & np.isfinite(voffs_tofit)

    model = modeling.powerlaws.PowerLaw1D(30, 20, 0.5)
    result = fitter(model, xoffs_au_tofit[ok], voffs_tofit[ok].to(u.km/u.s))

    model = modeling.powerlaws.PowerLaw1D(30, 20, 0.5)
    model.fixed['alpha'] = True
    result_alphapt5 = fitter(model, xoffs_au_tofit[ok], voffs_tofit[ok].to(u.km/u.s))

    mass_of_pt5 = (result_alphapt5.amplitude**2 / constants.G * result_alphapt5.x_0).to(u.M_sun)
    print("Best-fit mass with fixed alpha=0.5={0}".format(mass_of_pt5))

    model = modeling.powerlaws.PowerLaw1D(30, 20, 1)
    model.fixed['alpha'] = True
    result_alpha1 = fitter(model, xoffs_au_tofit[ok], voffs_tofit[ok].to(u.km/u.s))

    full_xarr = np.linspace(0, 100)*u.au
    pl.figure(1).clf()
    pl.xlabel("Offset from center (AU)")
    pl.ylabel("Offset from centroid velocity (km s$^{-1}$)")
    pl.plot(xoffs_au_tofit[ok & (xoffs_au>0)], voffs_tofit[ok & (xoffs_au>0)].to(u.km/u.s), 'o', color='r')
    pl.plot(xoffs_au_tofit[ok & (xoffs_au<0)], voffs_tofit[ok & (xoffs_au<0)].to(u.km/u.s), 's', color='b')
    pl.plot(full_xarr, result(full_xarr),
            label="$\\alpha={0:0.2f}$".format(result.alpha.value),
            color='k',
           )
    pl.plot(full_xarr,
            result_alpha1(full_xarr),
            label="$\\alpha=1$",
            linestyle='--',
            color='g',
           )
    pl.plot(full_xarr,
            result_alphapt5(full_xarr),
            label="$\\alpha=0.5$",
            linestyle='-.',
            color='m',
           )

    mass_lo = 13*u.M_sun
    mass_hi = 17*u.M_sun
    vel_lo = (((constants.G * mass_lo)/(full_xarr))**0.5).to(u.km/u.s)
    vel_hi = (((constants.G * mass_hi)/(full_xarr))**0.5).to(u.km/u.s)
    pl.fill_between(full_xarr.value,
                    vel_lo.value,
                    vel_hi.value,
                    zorder=-100, alpha=0.1,
                    color='m',
                    facecolor='m')


    pl.axis([0,90,12,30])
    pl.legend(loc='best')
    pl.savefig(paths.fpath("pv/bestfit_powerlaw_{0}".format(key)))

    pl.figure(2, figsize=(12,6)).clf()
    pl.subplot(1,3,1)
    pl.plot(xoffs_au_tofit[(xoffs_au>xmin)], np.abs(all_voff_results[key]['voffs_3'][(xoffs_au>xmin)].to(u.km/u.s)-vcen), '<', markeredgecolor='r', markerfacecolor='none', label='$3\sigma$')
    pl.plot(xoffs_au_tofit[(xoffs_au<-xmin)], np.abs(all_voff_results[key]['voffs_3'][(xoffs_au<-xmin)].to(u.km/u.s)-vcen), 'v', markeredgecolor='b', markerfacecolor='none')
    pl.subplot(1,3,2)
    pl.plot(xoffs_au_tofit[(xoffs_au>xmin)], np.abs(all_voff_results[key]['voffs_5'][(xoffs_au>xmin)].to(u.km/u.s)-vcen), 'o', color='r', markerfacecolor='r', label='$5\sigma$')
    pl.plot(xoffs_au_tofit[(xoffs_au<-xmin)], np.abs(all_voff_results[key]['voffs_5'][(xoffs_au<-xmin)].to(u.km/u.s)-vcen), 's', color='b', markerfacecolor='b')
    pl.subplot(1,3,3)
    pl.plot(xoffs_au_tofit[(xoffs_au>xmin)], np.abs(all_voff_results[key]['voffs_7'][(xoffs_au>xmin)].to(u.km/u.s)-vcen), '>', markeredgecolor='r', markerfacecolor='none', label='$7\sigma$')
    pl.plot(xoffs_au_tofit[(xoffs_au<-xmin)], np.abs(all_voff_results[key]['voffs_7'][(xoffs_au<-xmin)].to(u.km/u.s)-vcen), '^', markeredgecolor='b', markerfacecolor='none')

    for spn in (1,2,3):
        pl.subplot(1,3,spn)
        pl.plot(full_xarr,
                result_alphapt5(full_xarr),
                linestyle='-.',
                color='m',
               )

        pl.fill_between(full_xarr.value,
                        vel_lo.value,
                        vel_hi.value,
                        zorder=-100, alpha=0.1,
                        color='m',
                        facecolor='m')


        pl.xlabel("Offset from center (AU)")
        pl.ylabel("Offset from centroid velocity (km s$^{-1}$)")
        pl.axis([0,90,12,30])
        pl.legend(loc='best')

    pl.savefig(paths.fpath("pv/outerenvelope_velocity_thresholds_{0}".format(key)))

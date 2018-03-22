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
from line_point_offset import offset_to_point
import imp
import edge_on_ring_velocity_model
from constants import d_orion, vcen
import show_pv
import pylab as pl

imp.reload(edge_on_ring_velocity_model)
imp.reload(show_pv)

source = 'sourceI'
#diskycoord_list = pyregion.open(paths.rpath("{0}_disk_pvextract.reg"
#                                            .format(source)))[0].coord_list
#diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoord_list[jj],
#                                                     diskycoord_list[jj+1])
#                                    for jj in range(0,
#                                                    len(diskycoord_list),
#                                                    2)], unit=(u.deg,
#                                                               u.deg),
#                                   frame='fk5')
diskycoord_list = regions.read_ds9(paths.rpath("{0}_disk_pvextract.reg"
                                               .format(source)))
diskycoords = coordinates.SkyCoord([diskycoord_list[0].start,
                                    diskycoord_list[0].end])


#source = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
#                             unit=(u.hour, u.deg))
# fitted *disk* center from diskmodel
source = coordinates.SkyCoord(83.81048617*u.deg, -5.37516858*u.deg, frame='icrs')
source = coordinates.SkyCoord(regions.read_ds9(paths.rpath('sourceI_center.reg'))[0].center)
print(source.to_string('hmsdms'))
extraction_path = pvextractor.Path(diskycoords, width=0.01*u.arcsec)
origin = offset_to_point(source.ra.deg,
                         source.dec.deg,
                         extraction_path)*u.deg

all_voff_results = {}

for fn, vmin, vmax, savename, rms, radii in [#('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_robust0.5_diskpv.fits', -0.0005, 0.055,
                                             # 'H2O_kepler_SeifriedPlot.png', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust0.5_diskpv_0.01.fits', -0.0005, 0.048,
                                              'H2O_kepler_SeifriedPlot_0.01arcsec.pdf', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust0.5_diskpv_0.1.fits', -0.0005, 0.055,
                                              'H2O_kepler_SeifriedPlot_0.1arcsec.pdf', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.1.fits', -0.0005, 0.03,
                                              'H2O_kepler_SeifriedPlot_0.1arcsec_robust-2.pdf', 0.7*u.mJy, [10,100]),
                                             ('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.2.fits', -0.0005, 0.03,
                                              'H2O_kepler_SeifriedPlot_0.2arcsec_robust-2.pdf', 0.7*u.mJy, [10,100]),
                                             #('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_B6_robust-2_diskpv_0.01.fits', -0.0005, 0.048,
                                             # 'H2O_kepler_SeifriedPlot_0.01arcsec.pdf', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_29SiOv=0_5-4_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              '29SiOv0_5-4_kepler_SeifriedPlot.pdf', 2*u.mJy, [10,100]),
                                             ('pv/sourceI_SiS_12-11_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              'SiS_12-11_kepler_SeifriedPlot.pdf', 1*u.mJy, [30,200]),
                                             ('pv/sourceI_Unknown_1_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_1_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80]),
                                             ('pv/sourceI_Unknown_4_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_4_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80]),
                                             ('pv/sourceI_Unknown_5_B6_robust0.5_diskpv_0.01.fits', -0.005, 0.02,
                                              'Unknown_5_kepler_SeifriedPlot.pdf', 0.5*u.mJy, [30,80]),
                                             ('pv/sourceI_29SiOv=0_2-1_B3_robust-2_diskpv_0.01.fits', -0.05, 1,
                                              '29SiOv0_2-1_kepler_SeifriedPlot.pdf', 1*u.mJy, [10,100]),
                                             ('pv/sourceI_SiOv=1_5-4_B6_robust0.5_diskpv_0.01.fits', -0.01, 0.05,
                                              'SiOv1_5-4_kepler_SeifriedPlot.pdf', 2*u.mJy, [10,100]),
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
    vcen_ypix_left = int(np.round(ww.sub([2]).wcs_world2pix([(vcen-5*u.km/u.s).to(u.m/u.s).value], 0)[0][0]))
    vcen_ypix_right = int(np.round(ww.sub([2]).wcs_world2pix([(vcen+5*u.km/u.s).to(u.m/u.s).value], 0)[0][0]))
    xcen_xpix = int(np.round(ww.sub([1]).wcs_world2pix([0], 0)[0][0]))

    # rms ~ 1 mJy
    threshold = rms * 5

    def descend(xpix, ypix=vcen_ypix, direction=1):
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
    voffs = u.Quantity(ww.sub([2]).wcs_pix2world([descend(ii,
                                                          direction=-cdelt_sign if ii < xcen_xpix else cdelt_sign,
                                                          ypix=vcen_ypix_left if ii < xcen_xpix else vcen_ypix_right)
                                                  for ii in range(data.shape[1])], 0),
                       u.m/u.s).squeeze()
    all_voff_results[savename] = {'xoffs_as': xoffs_as, 'xoffs_au': xoffs_au, 'voffs': voffs}

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

    lines = show_pv.show_keplercurves(ax, origin, maxdist, vcen, masses=[15, ],
                                      linestyles=['-'], colors=['r'],
                                      radii={15: (radii, ('m','m'))},
                                      show_other_powerlaws=True
                                     )
    fig.savefig(paths.fpath('pv/{0}'.format(savename.replace(".pdf","_withalpha1.pdf"))),
                dpi=200,
                bbox_inches='tight')


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
    print("Best-fit mass with fixed alpha=0.5")

    model = modeling.powerlaws.PowerLaw1D(30, 20, 1)
    model.fixed['alpha'] = True
    result_alpha1 = fitter(model, xoffs_au_tofit[ok], voffs_tofit[ok].to(u.km/u.s))

    pl.figure(1).clf()
    pl.xlabel("Offset from center (AU)")
    pl.ylabel("Offset from centroid velocity (km s$^{-1}$)")
    pl.plot(xoffs_au_tofit[ok & (xoffs_au>0)], voffs_tofit[ok & (xoffs_au>0)].to(u.km/u.s), 'o', color='r')
    pl.plot(xoffs_au_tofit[ok & (xoffs_au<0)], voffs_tofit[ok & (xoffs_au<0)].to(u.km/u.s), 's', color='b')
    pl.plot(xoffs_au_tofit, result(xoffs_au_tofit),
            label="$\\alpha={0:0.2f}$".format(result.alpha.value),
            color='k',
           )
    pl.plot(xoffs_au_tofit,
            result_alpha1(xoffs_au_tofit),
            label="$\\alpha=1$",
            linestyle='--',
            color='g',
           )
    pl.plot(xoffs_au_tofit,
            result_alphapt5(xoffs_au_tofit),
            label="$\\alpha=0.5$",
            linestyle=':',
            color='m',
           )


    pl.axis([0,90,12,30])
    pl.legend(loc='best')
    pl.savefig(paths.fpath("pv/bestfit_powerlaw_{0}".format(key)))

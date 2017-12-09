import numpy as np
import paths
import pvextractor
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy import coordinates
import pyregion
from line_point_offset import offset_to_point
import imp
import show_pv
imp.reload(show_pv)

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


source = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                              unit=(u.hour, u.deg))
extraction_path = pvextractor.Path(diskycoords, width=0.1*u.arcsec)
origin = offset_to_point(source.ra.deg,
                         source.dec.deg,
                         extraction_path)*u.deg

for fn, vmin, vmax, savename, rms in [('pv/sourceI_H2Ov2=1_5(5,0)-6(4,3)_robust0.5_diskpv.fits', -0.0005, 0.055,
                                       'H2O_kepler_SeifriedPlot.png', 1*u.mJy),
                                      ('pv/sourceI_29SiOv=0_2-1_robust-2_diskpv.fits', -0.05, 1,
                                       'SiOv0_2-1_kepler_SeifriedPlot.png', 1*u.mJy),
                                      ('pv/sourceI_29SiOv=0_5-4_robust0.5_diskpv.fits', -0.01, 0.05,
                                       'SiOv0_5-4_kepler_SeifriedPlot.png', 2*u.mJy),
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

    vcen = 6.0*u.km/u.s
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

    fig,ax = show_pv.show_pv(data, ww,
                             origin, vrange=[-40,55], vcen=vcen,
                             imvmin=vmin, imvmax=vmax)

    ax.plot(xoffs_as, voffs, 'o-', transform=ax.get_transform('world'), markersize=3, markeredgecolor='b')
    maxdist=150*u.au
    show_pv.show_keplercurves(ax, origin, maxdist, vcen, masses=[19], linestyles=['-'])
    ax.set_aspect(2)


    fig.savefig(paths.fpath('pv/{0}'.format(savename)),
                dpi=200,
                bbox_inches='tight')

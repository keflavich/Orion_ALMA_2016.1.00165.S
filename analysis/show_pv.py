import numpy as np
import itertools
from astropy import constants
from astropy import units as u
from astropy import stats

import pylab as pl
from astropy import wcs
from mpl_toolkits.axes_grid1 import make_axes_locatable


from constants import d_orion

def show_pv(data, ww, origin, vrange, vcen, imvmin, imvmax, contour=False):
    

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

    fig = pl.figure(1, figsize=(8,8))
    fig.clf()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8],projection=ww)
    assert ww.wcs.cunit[1] == 'm/s' # this is BAD BAD BAD but necessary

    good_limits = (np.array((np.argmax(np.isfinite(data.max(axis=0))),
                             data.shape[1] -
                             np.argmax(np.isfinite(data.max(axis=0)[::-1])) - 1
                            ))
                   )
    leftmost_position = ww.wcs_pix2world(good_limits[0],
                                         vrange[0]*scalefactor,
                                         0)[0]*u.arcsec
    rightmost_position = ww.wcs_pix2world(good_limits[1],
                                          vrange[0]*scalefactor,
                                          0)[0]*u.arcsec
    assert rightmost_position > 0
    maxdist = ((rightmost_position)*d_orion).to(u.au, u.dimensionless_angles())
    assert maxdist > 0


    #imvmin,imvmax = (np.nanmin(data[plotted_slice]),
    #                 np.nanmax(data[plotted_slice]))
    im = ax.imshow(data, cmap='gray_r',
                   vmin=imvmin, vmax=imvmax*1.1,
                   interpolation='none')
    ax.set_xlabel("Offset [\"]")
    ax.set_ylabel("$V_{LSR}$ [km/s]")

    if contour:
        std_est = stats.mad_std(data, ignore_nan=True)
        print("Estimate of standard deviation: {0}".format(std_est))
        con = ax.contour(data, colors=[contour]*5,
                         levels=np.array([5,10,15,20,25])*std_est,
                         linewidths=0.75)
    else:
        con = None


    trans = ax.get_transform('world')
    length = (50*u.au / (d_orion)).to(u.deg, u.dimensionless_angles())
    endpoints_x = u.Quantity([0.1*u.arcsec, 0.1*u.arcsec+length]) + leftmost_position
    ax.plot(endpoints_x.to(u.arcsec),
            ([vrange[0]+2]*2*u.km/u.s).to(u.m/u.s),
            'r',
            transform=trans,
            zorder=100, linewidth=2)
    ax.text(endpoints_x.mean().value,
            (vrange[0]+3)*scalefactor,
            "50 au", color='r', transform=trans, ha='center')
    ax.plot(u.Quantity([leftmost_position, rightmost_position]).value,
            u.Quantity([vcen,vcen]).to(u.m/u.s).value, 'w--', transform=trans)


    ax.vlines(0, #origin.to(u.arcsec).value,
              (vrange[0]-5)*scalefactor,
              (vrange[1]+5)*scalefactor,
              color='w', linestyle='--', linewidth=2.0,
              alpha=0.6, transform=trans)

    ax.set_xlim(good_limits)

    ax.set_ylim(ww.wcs_world2pix(0,vrange[0]*scalefactor,0)[1],
                ww.wcs_world2pix(0,vrange[1]*scalefactor,0)[1])
    ax.set_xlim(good_limits)
    #ax2.set_xlim(good_limits)


    # ax.set_aspect(4)
    #ax.set_aspect(2*data.shape[1]/data.shape[0])
    aspect = 1*data.shape[1]/data.shape[0]
    ax.set_aspect(aspect)

    #ax.set_aspect('equal')

    ax.coords[1].set_format_unit(u.km/u.s)


    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # https://stackoverflow.com/a/18195921/814354
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=-0.05)

    #cax = fig.add_axes([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])

    #cb = pl.colorbar(im, cax=cax)
    #cb = pl.colorbar(mappable=im, ax=ax, pad=0.01, aspect=aspect*2)
    cb = pl.colorbar(mappable=im)
    pl.draw()
    bb = ax.bbox._bbox
    cb.ax.set_position([bb.x1+0.03, bb.y0, 0.05, bb.y1-bb.y0])

    ax.set_xlim(good_limits)

    return fig,ax,cb,con


def show_keplercurves(ax, origin, maxdist, vcen, masses=[15],
                      linestyles='-',
                      linewidths=None,
                      colors=['r'],
                      radii={15: ([30, 80], ['m', 'm'])},
                      yaxis_unit=u.m/u.s,
                      show_other_powerlaws=False,
                     ):

    trans = ax.get_transform('world')

    # Since we reset the CRVAL to be the origin, we don't need to
    # add origin as an offset any more
    # Rather than remove it, I'm just resetting it to zero...
    origin = 0*u.arcsec

    # overlay a Keplerian velocity curve
    positions = u.Quantity(np.linspace(0,maxdist,1000), u.au)

    if linewidths is None:
        linewidths = itertools.cycle([1])

    for mass,color,linestyle,linewidth in zip(masses,colors,linestyles,linewidths):
        # this is the 3d velocity, so assumes edge-on
        vel = (((constants.G * mass*u.M_sun)/(positions))**0.5).to(yaxis_unit)
        loc = (positions/(d_orion)).to(u.arcsec, u.dimensionless_angles())
        ax.plot((origin+loc).to(u.arcsec), (vcen+vel).to(yaxis_unit), linestyle=linestyle,
                color=color, linewidth=linewidth, alpha=1.0, transform=trans)
        #ax.plot((origin-loc).to(u.arcsec), (vcen+vel).to(yaxis_unit), 'b:', linewidth=1.0, alpha=1.0, transform=trans)
        #ax.plot((origin+loc).to(u.arcsec), (vcen-vel).to(yaxis_unit), 'b:', linewidth=1.0, alpha=1.0, transform=trans)
        ax.plot((origin-loc).to(u.arcsec), (vcen-vel).to(yaxis_unit), linestyle=linestyle,
                color=color, linewidth=linewidth, alpha=1.0, transform=trans)

        if show_other_powerlaws:
            vcurve = positions**-1 / (positions[200]**-1) * (vel[200])
            ax.plot((origin+loc).to(u.arcsec),
                    (vcen+vcurve).to(yaxis_unit),
                    color=color, linestyle='--',
                    transform=trans)
            ax.plot((origin-loc).to(u.arcsec),
                    (vcen-vcurve).to(yaxis_unit),
                    color=color, linestyle='--',
                    transform=trans)


    lines = []
    for mass in radii:
        for radius,color in zip(*radii[mass]):
            rad_as = (radius*u.au/d_orion).to(u.arcsec, u.dimensionless_angles())
            vel = (((constants.G * mass*u.M_sun)/(radius*u.au))**0.5).to(yaxis_unit)
            print("rad, vel: {0}, {1}".format(rad_as, vel))
            lines += ax.plot(u.Quantity([-rad_as, rad_as]),
                             (vcen+u.Quantity([-vel, vel])).to(yaxis_unit),
                             color=color, linestyle='--',
                             transform=trans)

    return lines

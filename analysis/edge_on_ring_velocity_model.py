import numpy as np
from astropy import constants
from astropy import units as u


def thindiskcurve(mass=20*u.M_sun, maxdist=100*u.au, yaxis_unit=u.km/u.s,
                  rmin=20*u.au, rmax=50*u.au,
                  npix=500):

    yy,xx = np.indices([npix,npix]) * maxdist * 2 / npix
    rr = ((yy-maxdist)**2 + (xx-maxdist)**2)**0.5
    theta = np.arctan2((yy-maxdist), (xx-maxdist))
    rr = u.Quantity(rr, u.au)

    v_abs = (((constants.G * u.Quantity(mass,u.M_sun))/(rr))**0.5).to(yaxis_unit)

    v_los = v_abs * np.cos(theta)

    mask1 = np.ones(rr.shape, dtype='float')
    mask1[~((rr > rmin) & (rr < rmax))] = np.nan

    x_plt = xx[0,:] - maxdist

    return (u.Quantity(x_plt, u.au),
            u.Quantity(np.nanmean((v_los * mask1), axis=0), yaxis_unit))

    #pl.clf()
    #x_plt = (xx[0,:] - npix/2)/npix/2 * maxdist/np.sqrt(2)
    #pl.plot(x_plt, v_los.mean(axis=0))
    #pl.plot(x_plt, np.nanmean((v_los * mask1), axis=0))

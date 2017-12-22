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

if __name__ == "__main__":

    import pylab as pl
    import paths


    fig1 = pl.figure(1)
    fig1.clf()
    ax = fig1.gca()

    xx,yy = thindiskcurve(mass=20*u.M_sun, rmin=20*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=20, 20 < r < 50", linestyle='-')

    xx,yy = thindiskcurve(mass=20*u.M_sun, rmin=30*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=20, 30 < r < 50", linestyle='-')

    xx,yy = thindiskcurve(mass=20*u.M_sun, rmin=30*u.au, rmax=80*u.au)
    ax.plot(xx, yy, label="M=20, 30 < r < 80", linestyle='-')

    xx,yy = thindiskcurve(mass=5*u.M_sun, rmin=20*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=5, 20 < r < 50", linestyle='--')

    xx,yy = thindiskcurve(mass=5*u.M_sun, rmin=30*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=5, 30 < r < 50", linestyle='--')

    xx,yy = thindiskcurve(mass=5*u.M_sun, rmin=30*u.au, rmax=80*u.au)
    ax.plot(xx, yy, label="M=5, 30 < r < 80", linestyle='--')


    ax.set_xlabel("Offset (AU)")
    ax.set_ylabel("$V_{obs}$ [km s$^{-1}$]")
    pl.legend(loc='best')
    pl.savefig(paths.fpath('velcentroid/radius_mass_demo.pdf'), bbox_inches='tight')

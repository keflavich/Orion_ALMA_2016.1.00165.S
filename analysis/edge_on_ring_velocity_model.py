import numpy as np
from astropy import constants
from astropy import units as u


def thindiskcurve(mass=20*u.M_sun, maxdist=100*u.au, yaxis_unit=u.km/u.s,
                  rmin=20*u.au, rmax=50*u.au,
                  vgrid=np.arange(-50,50,1)*u.km/u.s,
                  npix=1500, pvd=False):

    vstep = np.diff(vgrid).mean()

    yy,xx = np.indices([npix,npix]) * maxdist * 2 / npix
    rr = ((yy-maxdist)**2 + (xx-maxdist)**2)**0.5
    theta = np.arctan2((yy-maxdist), (xx-maxdist))
    rr = u.Quantity(rr, u.au)

    v_abs = (((constants.G * u.Quantity(mass,u.M_sun))/(rr))**0.5).to(yaxis_unit)

    v_los = v_abs * np.cos(theta)

    mask1 = np.ones(rr.shape, dtype='float')
    mask1[~((rr > rmin) & (rr < rmax))] = np.nan

    x_plt = xx - maxdist

    if pvd:
        xim = []
        for vv in vgrid:
            mask = (v_los > vv) & (v_los < vv+vstep) & np.isfinite(mask1)
            if np.any(mask):
                xim_ = ((mask).sum(axis=0))
                xim.append(xim_)
            else:
                xim.append(np.ones(v_los.shape[1])*np.nan)

        return x_plt[0,:], vgrid, np.array(xim)

    else:
        xpos = []
        for vv in vgrid:
            mask = (v_los > vv) & (v_los < vv+vstep) & np.isfinite(mask1)
            if np.any(mask):
                xxv = (v_los[mask] * x_plt[mask]).sum() / v_los[mask].sum()
                xpos.append(xxv)
            else:
                xpos.append(u.Quantity(np.nan, maxdist.unit))

        return u.Quantity(xpos, maxdist.unit), vgrid

    # this is the other version: mean velocity at each position instead
    # of mean position at each velocity
    return (u.Quantity(x_plt, u.au),
            u.Quantity(np.nanmean((v_los * mask1), axis=0), yaxis_unit))

if __name__ == "__main__":

    import pylab as pl
    import paths


    fig1 = pl.figure(1)
    fig1.clf()
    ax = fig1.gca()

    xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=20*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=15, 20 < r < 50", linestyle=':')

    xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=30*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=15, 30 < r < 50", linestyle=':')

    xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=30*u.au, rmax=80*u.au)
    ax.plot(xx, yy, label="M=15, 30 < r < 80", linestyle=':')

    xx,yy = thindiskcurve(mass=10*u.M_sun, rmin=20*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=10, 20 < r < 50", linestyle='--')

    xx,yy = thindiskcurve(mass=10*u.M_sun, rmin=30*u.au, rmax=50*u.au)
    ax.plot(xx, yy, label="M=10, 30 < r < 50", linestyle='--')

    xx,yy = thindiskcurve(mass=10*u.M_sun, rmin=30*u.au, rmax=80*u.au)
    ax.plot(xx, yy, label="M=10, 30 < r < 80", linestyle='--')

    #xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=30*u.au, rmax=70*u.au)
    #ax.plot(xx, yy, label="M=15, 30 < r < 70", linestyle=':')


    #xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=10*u.au, rmax=40*u.au)
    #ax.plot(xx, yy, label="M=15, 10 < r < 40", linestyle=':')


    ax.set_xlabel("Offset (AU)")
    ax.set_ylabel("$V_{obs}$ [km s$^{-1}$]")
    pl.legend(loc='best')
    pl.savefig(paths.fpath('velcentroid/radius_mass_demo_big.pdf'), bbox_inches='tight')


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

    #xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=30*u.au, rmax=70*u.au)
    #ax.plot(xx, yy, label="M=15, 30 < r < 70", linestyle=':')


    #xx,yy = thindiskcurve(mass=15*u.M_sun, rmin=10*u.au, rmax=40*u.au)
    #ax.plot(xx, yy, label="M=15, 10 < r < 40", linestyle=':')


    ax.set_xlabel("Offset (AU)")
    ax.set_ylabel("$V_{obs}$ [km s$^{-1}$]")
    pl.legend(loc='best')
    pl.savefig(paths.fpath('velcentroid/radius_mass_demo.pdf'), bbox_inches='tight')

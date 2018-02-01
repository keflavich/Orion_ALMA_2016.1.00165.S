import sys
import numpy as np
from astropy import constants
from astropy import units as u
import lmfit
from constants import vcen


def thindiskcurve(mass=20*u.M_sun, maxdist=100*u.au, yaxis_unit=u.km/u.s,
                  rmin=20*u.au, rmax=50*u.au,
                  vgrid=np.arange(-50,50,1)*u.km/u.s,
                  npix=1500, pvd=False):
    """
    Parameters
    ----------
    pvd : bool
        Return a position-velocity diagram?
    """

    mass = u.Quantity(mass, u.M_sun)
    rmin = u.Quantity(rmin, u.au)
    rmax = u.Quantity(rmax, u.au)

    if rmin >= rmax:
        return

    vsteps = np.diff(vgrid)
    vsteps = u.Quantity(np.concatenate([vsteps.value, [vsteps[-1].value]]), vgrid.unit)

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
        for vv,vstep in zip(vgrid, vsteps):
            mask = (v_los > vv-vstep/2) & (v_los < vv+vstep/2) & np.isfinite(mask1)
            if np.any(mask):
                xim_ = ((mask).sum(axis=0))
                xim.append(xim_)
            else:
                xim.append(np.ones(v_los.shape[1])*np.nan)

        return x_plt[0,:], vgrid, np.array(xim)

    else:
        xpos = []
        for vv,vstep in zip(vgrid, vsteps):
            mask = (v_los > vv-vstep/2) & (v_los < vv+vstep/2) & np.isfinite(mask1)
            if np.any(mask):
                if v_los[mask].sum() > 0.1*u.km/u.s:
                    xxv = (v_los[mask] * x_plt[mask]).sum() / v_los[mask].sum()
                else:
                    # around v=0, the velocity-weighted mean is infinite-ish
                    xxv = x_plt[mask].mean()
                xpos.append(xxv)
            else:
                xpos.append(u.Quantity(np.nan, maxdist.unit))

        return u.Quantity(xpos, maxdist.unit), vgrid

    # this is the other version: mean velocity at each position instead
    # of mean position at each velocity
    return (u.Quantity(x_plt, u.au),
            u.Quantity(np.nanmean((v_los * mask1), axis=0), yaxis_unit))

def thindiskcurve_residual(parameters, xsep, velo, error=None, **kwargs):

    # cheap progressbar
    sys.stdout.write('.')

    mass, rinner, delta, router, vcen = parameters.values()

    result = thindiskcurve(mass=mass,
                           rmin=rinner,
                           rmax=router,
                           vgrid=velo-u.Quantity(vcen.value, u.km/u.s),
                           **kwargs)
    if result is None:
        return np.ones_like(xsep) * 1e10

    model_xsep, model_v = result

    resid = (xsep-model_xsep).value

    if error is None:
        error = np.ones_like(resid)

    #error[np.isnan(resid)] = 1e10
    # we want significant, but hopefully not dominant, residuals when the model
    # predicts no data but there are data
    resid[np.isnan(resid)] = 15

    return resid/error

def thindiskcurve_fitter(xsep, velo, error=None, mguess=20*u.M_sun,
                         rinner=20*u.au, router=50*u.au):

    parameters = lmfit.Parameters()
    parameters.add('mass', value=u.Quantity(mguess, u.M_sun).value, min=10, max=25)
    parameters.add('rinner', value=u.Quantity(rinner, u.au).value, min=3, max=50)
    parameters.add('delta', value=20, min=10, max=50)
    parameters.add('router', value=u.Quantity(router, u.au).value, min=20, max=100,
                   expr='rinner+delta')
    parameters.add('vcen', value=vcen.value, min=3.5, max=7.5)
    result = lmfit.minimize(thindiskcurve_residual, parameters, epsfcn=0.005,
                            kws={'xsep': u.Quantity(xsep, u.au),
                                 'velo': u.Quantity(velo, u.km/u.s),
                                 'error': error})

    result.params.pretty_print()

    return result

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

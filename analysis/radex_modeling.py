import pyradex
import pylab as pl
import numpy as np
import paths
import dust_emissivity
from astropy import units as u
from astropy.table import Table, Column
from astropy.modeling import models, fitting
import lmfit

rr = pyradex.Radex(species='nacl', temperature=1000, density={'H2': 1e8}, column=4e13)
rslt = rr()

tbl = Table.read(paths.tpath('line_fits.txt'), format='ascii.fixed_width')
tbl.add_column(name='RADEX_row', col=Column(np.zeros(len(tbl), dtype='int')))

for row in tbl:
    upper = ("{0:<6}".format(f"{row['v']}_{row['J$_u$']}")).encode()
    lower = ("{0:<6}".format(f"{row['v']}_{row['J$_l$']}")).encode()
    match = (rslt['upperlevel'] == upper) & (rslt['lowerlevel'] == lower)
    if row['v'] <= 3:
        assert match.sum() == 1
    if match.sum() == 1:
        row['RADEX_row'] = np.argmax(match)
    else:
        row['RADEX_row'] = -999

fittbl_cache = {}

def get_fittbl(species):

    if species not in fittbl_cache:

        halogen = species.lower().split('cl')[0]
        fittbl = tbl[(tbl['RADEX_row'] >= 0) &
                     (tbl['Fitted Amplitude error K'] < tbl['Fitted Amplitude K']) &
                     (np.array([halogen in row['Line Name'].lower() for row in tbl]))
                    ]

        fittbl_cache[species] = fittbl

    return fittbl_cache[species]

def resid(rdx, species='nacl'):

    tbl = get_fittbl(species)

    inds = tbl['RADEX_row']

    # need non-background-subtracted brightness
    predicted_brightness = rdx['Tex'][inds] * (1-np.exp(-rdx['tau'][inds]))

    obs = tbl['Fitted Amplitude K']
    err = tbl['Fitted Amplitude error K']

    return (((obs-predicted_brightness)/err))

def chi2(rdx, species='nacl'):
    return (resid(rdx, species=species)**2).sum()

def GOF_plot(rdx, species='nacl'):

    tbl = get_fittbl(species)

    inds = tbl['RADEX_row']
    predicted_brightness = rdx['Tex'][inds] * (1-np.exp(-rdx['tau'][inds]))

    obs = tbl['Fitted Amplitude K']
    err = tbl['Fitted Amplitude error K']

    energy = rdx['upperstateenergy'][inds]

    pl.plot(energy, predicted_brightness, 's')
    pl.errorbar(energy, obs, yerr=err, linestyle='none', marker='.')


freq = rr.frequency
wl = freq.to(u.um, u.spectral())
rovib_range = ((wl>25*u.um) & (wl<45*u.um))

def bgfunc(freq, disktem=150*u.K, diskdilution=1/2., startem=4000*u.K,
           stardilution=(100*u.R_sun)**2/(30*u.au)**2, cmbtem=2.73*u.K,
           cmbdilution=1.0):
    return (dust_emissivity.blackbody.blackbody(nu=freq, temperature=disktem) * diskdilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=startem) * stardilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=cmbtem) * cmbdilution)

def bg_with_artificial(freq, rovib_range=rovib_range, gcen=29*u.um,
                       gwidth=1.0*u.um,
                       gamp=5e-8*rr.background_brightness.unit, **kwargs):
    wl = freq.to(u.um, u.spectral())
    artificial = bgfunc(freq, **kwargs)
    artificial[rovib_range] += gamp * np.exp(-(wl[rovib_range]-gcen)**2/(2*(gwidth)**2))
    return artificial

def fitmod(gcen=29*u.um, gwidth=1*u.um, gamp=5e-8*rr.background_brightness.unit,
           density=1e8*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K,
           disktem=150*u.K, diskdilution=1/2., startem=4000*u.K,
           stardilution=(100*u.R_sun)**2/(30*u.au)**2,
           plot=False, species='nacl'
          ):

    gcen = u.Quantity(gcen, u.um)
    gwidth = u.Quantity(gwidth, u.um)
    gamp = u.Quantity(gamp, rr.background_brightness.unit)
    density = u.Quantity(density, u.cm**-3)
    column = u.Quantity(column, u.cm**-2)
    temperature = u.Quantity(temperature, u.K)
    disktem = u.Quantity(disktem, u.K)
    startem = u.Quantity(startem, u.K)
    stardilution = u.Quantity(stardilution, u.dimensionless_unscaled)
    diskdilution = u.Quantity(diskdilution, u.dimensionless_unscaled)

    rr.background_brightness = bg_with_artificial(freq, gcen=gcen,
                                                  gwidth=gwidth, gamp=gamp,
                                                  diskdilution=diskdilution,
                                                  stardilution=stardilution,
                                                  disktem=disktem,
                                                  startem=startem,
                                                 )
    rslt2 = (rr(density={'H2':density}, column=column, temperature=temperature,
                species=species,
                tbg=None))

    if plot:
        GOF_plot(rslt2, species=species)

    return resid(rslt2, species=species)

def lmfitmod(kwargs):
    return fitmod(**{key: par.value for key, par in kwargs.items()})

def iter_cb(params, iter, resid, *args, **kws):
    formatted_pars = " ".join([f"{key}:{value.value:0.5g}" for key,value in
                               params.items()])
    print(f"iter={iter:04d} params={formatted_pars} chi2={(resid**2).sum():10.2f}")


# defaults
pars = lmfit.Parameters()
pars.add(name='gcen', value=29, min=27, max=31)
pars.add(name='gwidth', value=1, min=0.5, max=1.5)
pars.add(name='gamp', value=5e-8, min=1e-12, max=1e-5)
pars.add(name='density', value=1e8, min=1e6, max=1e10)
pars.add(name='column', value=1e14, min=1e12, max=1e17)
pars.add(name='temperature', value=100, min=50, max=2000)

# previous best fit
pars = lmfit.Parameters()
pars.add(name='gcen', value=28.981, min=27.0, max=31.0)
pars.add(name='gwidth', value=0.309, min=0.25, max=1.5)
pars.add(name='gamp', value=1.62e-5, min=1e-12, max=1e-4)
pars.add(name='density', value=8.06e7, min=1e6, max=1e10)
pars.add(name='column', value=9.11e13, min=1e12, max=1e17)
pars.add(name='temperature', value=904, min=50, max=2000)

def fit():
    return lmfit.minimize(lmfitmod, pars, iter_cb=iter_cb)

noIRpars = lmfit.Parameters()
#noIRpars.add(name='gcen', value=28.9, min=27.0, max=31.0, vary=False)
#noIRpars.add(name='gwidth', value=0.0, min=0, vary=False)
#noIRpars.add(name='gamp', value=0, min=0, vary=False)
noIRpars.add(name='density', value=1e8, min=1e6, max=1e18)
noIRpars.add(name='column', value=1.16e14, min=1e12, max=1e17)
noIRpars.add(name='temperature', value=800, min=50, max=5000)
def fit_no_ir_line():
    return lmfit.minimize(lmfitmod, noIRpars, iter_cb=iter_cb)

starnoIRpars = lmfit.Parameters()
#starnoIRpars.add(name='gcen', value=28.9, min=27.0, max=31.0, vary=False)
#starnoIRpars.add(name='gwidth', value=0.0, min=0, vary=False)
#starnoIRpars.add(name='gamp', value=0, min=0, vary=False)
starnoIRpars.add(name='density', value=1.4e7, min=1e6, max=1e18)
starnoIRpars.add(name='column', value=6.1e13, min=1e12, max=1e17)
starnoIRpars.add(name='temperature', value=950, min=50, max=5000)
starnoIRpars.add(name='startem', value=4000.0, min=1000, max=10000)
starnoIRpars.add(name='disktem', value=150.0, min=50, max=400)
starnoIRpars.add(name='stardilution', value=2e-4, min=2e-9, max=2e-3)
starnoIRpars.add(name='diskdilution', value=0.5, min=0.1, max=1)
def fit_star_no_ir_line():
    return lmfit.minimize(lmfitmod, starnoIRpars, iter_cb=iter_cb)

import pyradex
import pylab as pl
import numpy as np
import paths
import dust_emissivity
from astropy import units as u
from astropy.table import Table
from kcl_rotation_diagram import fit_tex
from astroquery.vamdc import Vamdc
from radex_modeling import chi2, GOF_plot, fitmod

rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
rslt = rr()

obs = (((85.5 < rslt['frequency']) & (89.5 > rslt['frequency'])) |
       ((97.5 < rslt['frequency']) & (101.5 > rslt['frequency'])) |
       ((229.0 < rslt['frequency']) & (233.7 > rslt['frequency'])) |
       ((214.25 < rslt['frequency']) & (218.9 > rslt['frequency'])) |
       ((344.1 < rslt['frequency']) & (348 > rslt['frequency'])) |
       ((332.1 < rslt['frequency']) & (335.8 > rslt['frequency'])))


freq = rr.frequency
wl = freq.to(u.um, u.spectral())
rovib_range = ((wl>25*u.um) & (wl<45*u.um))

def bgfunc(freq, disktem=150*u.K, diskdilution=1/2., startem=4000*u.K, stardilution=(100*u.R_sun)**2/(30*u.au)**2,
           cmbtem=2.73*u.K, cmbdilution=1.0):
    return (dust_emissivity.blackbody.blackbody(nu=freq, temperature=disktem) * diskdilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=startem) * stardilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=cmbtem) * cmbdilution)

def bg_with_artificial(freq, rovib_range=rovib_range, gcen=29*u.um,
                       gwidth=1.0*u.um,
                       gamp=5e-8*rr.background_brightness.unit, **kwargs):
    wl = freq.to(u.um, u.spectral())
    artificial = bgfunc(freq)
    artificial[rovib_range] += gamp * np.exp(-(wl[rovib_range]-gcen)**2/(2*(gwidth)**2))
    return artificial

rr.background_brightness = bg_with_artificial(freq)


plotwl = np.logspace(0, 5, 1000)*u.um
plotfreq = plotwl.to(u.GHz, u.spectral())
rovib_range_plot = (plotwl>25*u.um) & (plotwl<45*u.um)
plotbg = bg_with_artificial(plotfreq, rovib_range_plot)

rslt = (rr(density=1e8*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None))
print("wacky background-illuminated chi^2 = {0}".format(chi2(rslt)))
vone = np.array([row['upperlevel'][0] == '1' for row in rslt], dtype='bool')
vzero = np.array([row['upperlevel'][0] == '0' for row in rslt], dtype='bool')
vtwo = np.array([row['upperlevel'][0] == '2' for row in rslt], dtype='bool')
Jeight = np.array([row['upperlevel'][2] == '8' for row in rslt], dtype='bool')
pl.figure(1).clf()
pl.subplot(2,1,1)
pl.loglog(plotwl, plotbg, linestyle='-')
pl.loglog(wl, rr.background_brightness, marker='.', linestyle='none')
pl.xlabel("Wavelength [$\mu$m]")
pl.ylabel("Background Brightness\n[{}]".format(rr.background_brightness.unit.to_string(format='latex')))
pl.subplot(2,1,2)

# introduce an arbitrary scaling value to obtain 'reasonable' column densities
scaling = 1e12 * u.cm**-2

pl.plot(rslt['upperstateenergy'].data, np.log10(rslt['upperlevelpop'].data * scaling.value), '.', label=None, alpha=0.15)
#pl.plot(rslt[vzero]['upperstateenergy'].data, np.log10(rslt[vzero]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[vone]['upperstateenergy'].data, np.log10(rslt[vone]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[vzero & obs]['upperstateenergy'].data, np.log10(rslt[vzero & obs]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[vone & obs]['upperstateenergy'].data, np.log10(rslt[vone & obs]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[vtwo]['upperstateenergy'].data, np.log10(rslt[vtwo]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[vtwo & obs]['upperstateenergy'].data, np.log10(rslt[vtwo & obs]['upperlevelpop'].data), 'o', label=None)
#pl.plot(rslt[Jeight]['upperstateenergy'].data, np.log10(rslt[Jeight]['upperlevelpop'].data), 's', label=None)


nacl = Vamdc.query_molecule('NaCl$')

fit_tex(u.Quantity(rslt[vone & obs]['upperstateenergy'], u.K), rslt[vone &
                                                                    obs]['upperlevelpop']*scaling,
        molecule=nacl,plot=True, min_nupper=1e-50)
fit_tex(u.Quantity(rslt[Jeight]['upperstateenergy'], u.K),
        rslt[Jeight]['upperlevelpop']*scaling, molecule=nacl,plot=True,
        min_nupper=1e-50, color='k')
pl.legend(loc='upper right')

pl.xlabel("E$_U$ [K]")
pl.ylabel("log upper state population")
pl.tight_layout()
pl.savefig(paths.fpath('simulated_populations_with_wacky_radiation_field.pdf'))

# the referee asked about tau
# blue = all
# orange = observed
pl.figure(2).clf()
pl.plot(rslt['upperstateenergy'].data, rslt['tau'], '.')
pl.plot(rslt['upperstateenergy'].data[obs], rslt['tau'][obs], '.')
pl.xlabel("Upper state energy")
pl.ylabel("Optical Depth")

# What about the observable T_B?
# (note that RADEX is using the background-subtracted brightness in the T_B column, which isn't what we want)
pl.figure(3).clf()
T_B = (rr.source_brightness*u.sr).to(u.K, u.brightness_temperature(1*u.sr, rr.frequency))
pl.plot(rslt['upperstateenergy'].data, T_B, '.')
pl.plot(rslt['upperstateenergy'].data[obs], T_B[obs], '.')
pl.xlabel("Upper state energy")
pl.ylabel("Brightness Temperature")
pl.ylim(0,100)




pl.figure(4).clf()
GOF_plot(rslt)

pl.figure(5).clf()
rr.background_brightness = bg_with_artificial(freq, gcen=29*u.um, gwidth=1*u.um,
                                              gamp=5e-8*rr.background_brightness.unit,)
rslt2 = (rr(density=1e8*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=None))
GOF_plot(rslt2)

pl.figure(6).clf()
fitmod(plot=True)

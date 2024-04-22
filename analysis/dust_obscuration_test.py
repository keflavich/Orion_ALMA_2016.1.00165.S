"""
The referee suggested that dust attenuation in the upper disk might be
responsible for the lower intensities of the high-J lines, i.e., the
apparently cool rotational temperatures.  This probably can't be a complete
explanation, since there still needs to be very high excitation to get
the v=n states occupied, but we can evaluate this idea.


This wasn't the most helpful test, though, because it doesn't go back through
the step of inferring the level populations afterward.
"""
import pyradex
import pylab as pl
import numpy as np
import paths
import dust_emissivity
from astropy import units as u
from astropy.table import Table, Column
from kcl_rotation_diagram import fit_tex
#from astroquery.vamdc import Vamdc
from dust_emissivity.dust import kappa
from radex_modeling import chi2, GOF_plot
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters

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

def bgfunc(freq, disktem=150*u.K, diskdilution=1/2.,
           startem=4000*u.K,
           stardilution=(100*u.R_sun)**2/(30*u.au)**2,
           cmbtem=2.73*u.K, cmbdilution=1.0):
    return (dust_emissivity.blackbody.blackbody(nu=freq, temperature=disktem) * diskdilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=startem) * stardilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=cmbtem) * cmbdilution)

# no artificial background, just a reasonable one
rr.background_brightness = bgfunc(freq)


plotwl = np.logspace(0, 5, 1000)*u.um
plotfreq = plotwl.to(u.GHz, u.spectral())
rovib_range_plot = (plotwl>25*u.um) & (plotwl<45*u.um)
plotbg = bgfunc(plotfreq)

rslt = (rr(density=1e8*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None))
print("stellar background-illuminated chi^2 = {0}".format(chi2(rslt)))
vone = np.array([row['upperlevel'][0] == '1' for row in rslt], dtype='bool')
vzero = np.array([row['upperlevel'][0] == '0' for row in rslt], dtype='bool')
vtwo = np.array([row['upperlevel'][0] == '2' for row in rslt], dtype='bool')
Jeight = np.array([row['upperlevel'][2] == '8' for row in rslt], dtype='bool')
pl.figure(1)
pl.clf()
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


frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('NaCl, v=0-15', catalog='CDMS', fmin=85*u.GHz, fmax=360*u.GHz)

fit_tex(u.Quantity(rslt[vone & obs]['upperstateenergy'], u.K), rslt[vone &
                                                                    obs]['upperlevelpop']*scaling,
        partition_func=partfunc, plot=True, min_nupper=1e-50)
fit_tex(u.Quantity(rslt[Jeight]['upperstateenergy'], u.K),
        rslt[Jeight]['upperlevelpop']*scaling, partition_func=partfunc,
        plot=True,
        min_nupper=1e-50, color='k')
pl.legend(loc='upper right')

pl.xlabel("E$_U$ [K]")
pl.ylabel("log upper state population")
pl.tight_layout()
pl.savefig(paths.fpath('simulated_populations_with_stellar_radiation_field.pdf'))



# the referee asked about tau
# blue = all
# orange = observed
pl.figure(2).clf()
pl.plot(rslt['upperstateenergy'].data, rslt['tau'], '.')
pl.plot(rslt['upperstateenergy'].data[obs], rslt['tau'][obs], '.')
pl.xlabel("Upper state energy")
pl.ylabel("Optical Depth")

pl.figure(3).clf()

T_B = (rr.source_brightness*u.sr).to(u.K, u.brightness_temperature(rr.frequency, 1*u.sr, ))

# what if we attenuate the emission either by increasing the background or
# decreasing the observed emission?

T_B = T_B - 100 * u.K * (1-np.exp(-kappa(rr.frequency, beta=1)*10*u.g/u.cm**2))

pl.plot(rslt['upperstateenergy'].data, T_B, '.')
pl.plot(rslt['upperstateenergy'].data[obs], T_B[obs], '.')
pl.xlabel("Upper state energy")
pl.ylabel("Brightness Temperature")
pl.ylim(0,100)


pl.figure(4).clf()
GOF_plot(rslt)

pl.figure(5).clf()
rslt2 = rr(density=1e8*u.cm**-3, column=1e19*u.cm**-2, temperature=90*u.K, tbg=None)
GOF_plot(rslt2)

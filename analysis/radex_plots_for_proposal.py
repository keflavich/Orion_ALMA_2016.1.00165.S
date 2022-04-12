import pyradex
import pylab as pl
import numpy as np
import paths
import dust_emissivity
from astropy import units as u
from astropy.table import Table
from radex_modeling import chi2

rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
rslt = rr()


v0_76 = (rslt['upperlevel'] == '0_7   ') & (rslt['lowerlevel'] == '0_6   ')
v1_76 = (rslt['upperlevel'] == '1_7   ') & (rslt['lowerlevel'] == '1_6   ')
v2_76 = (rslt['upperlevel'] == '2_7   ') & (rslt['lowerlevel'] == '2_6   ')
v3_76 = (rslt['upperlevel'] == '3_7   ') & (rslt['lowerlevel'] == '3_6   ')
v10 = (rslt['upperlevel'] == '1_7   ') & (rslt['lowerlevel'] == '0_6   ')
v21 = (rslt['upperlevel'] == '2_7   ') & (rslt['lowerlevel'] == '1_6   ')
v32 = (rslt['upperlevel'] == '3_7   ') & (rslt['lowerlevel'] == '2_6   ')


v0_1817 = (rslt['upperlevel'] == '0_18  ') & (rslt['lowerlevel'] == '0_17  ')
v1_1817 = (rslt['upperlevel'] == '1_18  ') & (rslt['lowerlevel'] == '1_17  ')
v2_1817 = (rslt['upperlevel'] == '2_18  ') & (rslt['lowerlevel'] == '2_17  ')
v3_1817 = (rslt['upperlevel'] == '3_18  ') & (rslt['lowerlevel'] == '3_17  ')

obs = (((89.5 < rslt['frequency']) & (93.5 > rslt['frequency'])) |
       ((101.5 < rslt['frequency']) & (105.25 > rslt['frequency'])) |
       ((153 < rslt['frequency']) & (157 > rslt['frequency'])) |
       ((141 < rslt['frequency']) & (145 > rslt['frequency'])) |
       ((244 < rslt['frequency']) & (248 > rslt['frequency'])) |
       ((257.5 < rslt['frequency']) & (261.5 > rslt['frequency'])) |
       ((348.8 < rslt['frequency']) & (352.3 > rslt['frequency'])) |
       ((361 < rslt['frequency']) & (364 > rslt['frequency'])) |
       ((479.1 < rslt['frequency']) & (482.9 > rslt['frequency'])) |
       ((491.2 < rslt['frequency']) & (493 > rslt['frequency']))
      )


# play with different backgrounds

freq = rr.frequency
wl = rr.frequency.to(u.um, u.spectral())
bb100_half_plus_cmb = dust_emissivity.blackbody.blackbody(nu=freq, temperature=100*u.K)/2. + dust_emissivity.blackbody.blackbody(nu=freq, temperature=2.73*u.K)

rr.background_brightness = bb100_half_plus_cmb

bb4000smallff_100_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=100*u.K)/2. +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=2.73*u.K))

rr.background_brightness = bb4000smallff_100_half_plus_cmb

bb4000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=200*u.K)/2. +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=2.73*u.K))
rr.background_brightness = bb4000smallff_200_half_plus_cmb

mbb1000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=200*u.K)/2. +
                                    dust_emissivity.blackbody.modified_blackbody(nu=freq,
                                                                                 temperature=1000*u.K)*(10*u.au)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=2.73*u.K))
rr.background_brightness = mbb1000smallff_200_half_plus_cmb

rovib_range = ((wl>25*u.um) & (wl<45*u.um))
stepfunc = 1e-7*u.erg/u.s/u.cm**2/u.Hz/u.sr*rovib_range
artificial = bb4000smallff_200_half_plus_cmb + stepfunc
rr.background_brightness = artificial
print(rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))


rr.background_brightness = mbb1000smallff_200_half_plus_cmb
rr.background_brightness = artificial
rr.background_brightness = bb100_half_plus_cmb

artificial_two = bb100_half_plus_cmb.copy()/4
artificial_two[rovib_range] = 1e-10*rr.background_brightness.unit * (np.arange(1, rovib_range.sum()+1)/rovib_range.sum()*4 + 1)[::-1]
rr.background_brightness = artificial_two

def bgfunc(freq, disktem=100*u.K, diskdilution=1/5., startem=4000*u.K, stardilution=(100*u.R_sun)**2/(30*u.au)**2, 
           cmbtem=2.73*u.K, cmbdilution=1.0):
    return (dust_emissivity.blackbody.blackbody(nu=freq, temperature=disktem) * diskdilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=startem) * stardilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=cmbtem) * cmbdilution)
artificial_three = bgfunc(freq)
artificial_three[rovib_range] += 1e-9*rr.background_brightness.unit * np.exp(-(wl[rovib_range]-35*u.um)**2/(2*(2.5*u.um)**2))
rr.background_brightness = artificial_three


plotwl = np.logspace(0, 5, 1000)*u.um
plotfreq = plotwl.to(u.GHz, u.spectral())
plotbg = bgfunc(plotfreq)
rovib_range_plot = (plotwl>25*u.um) & (plotwl<45*u.um)
plotbg[rovib_range_plot] += 1e-9*rr.background_brightness.unit * np.exp(-(plotwl[rovib_range_plot]-35*u.um)**2/(2*(2.5*u.um)**2))

rslt = (rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None))
print(chi2(rslt))
vone = np.array([row['upperlevel'][0] == '1' for row in rslt], dtype='bool')
vzero = np.array([row['upperlevel'][0] == '0' for row in rslt], dtype='bool')
vtwo = np.array([row['upperlevel'][0] == '2' for row in rslt], dtype='bool')
Jeight = np.array([row['upperlevel'][2] == '8' for row in rslt], dtype='bool')
pl.clf()
pl.subplot(2,1,1)
pl.loglog(plotwl, plotbg, linestyle='-')
pl.loglog(wl, rr.background_brightness, marker='.', linestyle='none')
pl.xlabel("Wavelength [$\mu$m]")
pl.ylabel("Background Brightness\n[{}]".format(rr.background_brightness.unit.to_string()))
pl.subplot(2,1,2)
pl.semilogy(rslt[vzero]['upperstateenergy'], rslt[vzero]['upperlevelpop'], 'o')
pl.semilogy(rslt[vone]['upperstateenergy'], rslt[vone]['upperlevelpop'], 'o')
pl.semilogy(rslt[vzero & obs]['upperstateenergy'], rslt[vzero & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[vone & obs]['upperstateenergy'], rslt[vone & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[vtwo]['upperstateenergy'], rslt[vtwo]['upperlevelpop'], 'o')
pl.semilogy(rslt[vtwo & obs]['upperstateenergy'], rslt[vtwo & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[Jeight]['upperstateenergy'], rslt[Jeight]['upperlevelpop'], 's')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Upper state population")
pl.tight_layout()
# MOVED to dust_obscuration pl.savefig(paths.fpath('simulated_populations_with_wacky_radiation_field.pdf'))


print(rr(density=1e4*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=1000*u.K)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
print(chi2(rr.get_table()))



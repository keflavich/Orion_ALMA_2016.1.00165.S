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

B3 = (((89.5 < rslt['frequency']) & (93.5 > rslt['frequency'])) |
       ((101.5 < rslt['frequency']) & (105.25 > rslt['frequency'])) )
B4 = (((153 < rslt['frequency']) & (157 > rslt['frequency'])) |
       ((141 < rslt['frequency']) & (145 > rslt['frequency'])) )
B6 = (((244 < rslt['frequency']) & (248 > rslt['frequency'])) |
       ((257.5 < rslt['frequency']) & (261.5 > rslt['frequency'])) )
B7 = (((348.8 < rslt['frequency']) & (352.3 > rslt['frequency'])) |
       ((361 < rslt['frequency']) & (364 > rslt['frequency'])) )
B8 = (
       ((479.1 < rslt['frequency']) & (482.9 > rslt['frequency'])) |
       ((491.2 < rslt['frequency']) & (493 > rslt['frequency']))
      )
obs= B3|B4|B6|B7|B8


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

#rr.background_brightness = bb4000smallff_100_half_plus_cmb
#
bb4000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=200*u.K)/2. +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=2.73*u.K))
#rr.background_brightness = bb4000smallff_200_half_plus_cmb
#
mbb1000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=200*u.K)/2. +
                                    dust_emissivity.blackbody.modified_blackbody(nu=freq,
                                                                                 temperature=1000*u.K)*(10*u.au)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=2.73*u.K))
#rr.background_brightness = mbb1000smallff_200_half_plus_cmb
#
rovib_range = ((wl>25*u.um) & (wl<45*u.um))
stepfunc = 1e-7*u.erg/u.s/u.cm**2/u.Hz/u.sr*rovib_range
artificial = bb4000smallff_200_half_plus_cmb + stepfunc
#rr.background_brightness = artificial
#print(rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
#print(chi2(rr.get_table()))
#
#
#rr.background_brightness = mbb1000smallff_200_half_plus_cmb
#rr.background_brightness = artificial
#rr.background_brightness = bb100_half_plus_cmb
#
artificial_two = bb100_half_plus_cmb.copy()/4
artificial_two[rovib_range] = 1e-10*rr.background_brightness.unit * (np.arange(1, rovib_range.sum()+1)/rovib_range.sum()*4 + 1)[::-1]
#rr.background_brightness = artificial_two

def bgfunc(freq, disktem=100*u.K, diskdilution=1/5., startem=4000*u.K, stardilution=(100*u.R_sun)**2/(30*u.au)**2, 
           cmbtem=2.73*u.K, cmbdilution=1.0):
    return (dust_emissivity.blackbody.blackbody(nu=freq, temperature=disktem) * diskdilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=startem) * stardilution +
            dust_emissivity.blackbody.blackbody(nu=freq, temperature=cmbtem) * cmbdilution)
artificial_three = bgfunc(freq)
artificial_three[rovib_range] += 1e-9*rr.background_brightness.unit * np.exp(-(wl[rovib_range]-35*u.um)**2/(2*(2.5*u.um)**2))
rr.background_brightness = artificial_three

artificial_29um = bgfunc(freq)
artificial_29um[rovib_range] += 1e-9*rr.background_brightness.unit * np.exp(-(wl[rovib_range]-29*u.um)**2/(2*(2.5*u.um)**2))


plotwl = np.logspace(0, 5, 1000)*u.um
plotfreq = plotwl.to(u.GHz, u.spectral())
plotbg = bgfunc(plotfreq)
rovib_range_plot = (plotwl>25*u.um) & (plotwl<45*u.um)
plotbg[rovib_range_plot] += 1e-9*rr.background_brightness.unit * np.exp(-(plotwl[rovib_range_plot]-35*u.um)**2/(2*(2.5*u.um)**2))

rslt = (rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=100*u.K, tbg=None))
#print(chi2(rslt))
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
pl.semilogy(rslt[vzero]['upperstateenergy'], rslt[vzero]['upperlevelpop'], ',')
pl.semilogy(rslt[vone]['upperstateenergy'], rslt[vone]['upperlevelpop'], ',')
pl.semilogy(rslt[vzero & obs]['upperstateenergy'], rslt[vzero & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[vone & obs]['upperstateenergy'], rslt[vone & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[vtwo]['upperstateenergy'], rslt[vtwo]['upperlevelpop'], ',')
pl.semilogy(rslt[vtwo & obs]['upperstateenergy'], rslt[vtwo & obs]['upperlevelpop'], 'o')
pl.semilogy(rslt[Jeight]['upperstateenergy'], rslt[Jeight]['upperlevelpop'], 's')
pl.xlabel("E$_U$ [K]")
pl.ylabel("Upper state population")
pl.tight_layout()
# MOVED to dust_obscuration pl.savefig(paths.fpath('simulated_populations_with_wacky_radiation_field.pdf'))


#print(rr(density=1e4*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=1000*u.K)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
#print(chi2(rr.get_table()))

def makeplot(rslt):
    ax1 = pl.subplot(2,1,1)
    L0, = ax1.semilogy(rslt[vzero]['upperstateenergy'], rslt[vzero]['T_B'], ',')
    L1, = ax1.semilogy(rslt[vone]['upperstateenergy'], rslt[vone]['T_B'], ',')
    L2, =ax1.semilogy(rslt[vtwo]['upperstateenergy'], rslt[vtwo]['T_B'], ',')
    ax1.semilogy(rslt[vzero & obs]['upperstateenergy'], rslt[vzero & obs]['T_B'], 'o', color=L0.get_color(), label='v=0')
    ax1.semilogy(rslt[vone & obs]['upperstateenergy'], rslt[vone & obs]['T_B'], 'o', color=L1.get_color(), label='v=1')
    ax1.semilogy(rslt[vtwo & obs]['upperstateenergy'], rslt[vtwo & obs]['T_B'], 'o', color=L2.get_color(), label='v=2')
    ax1.semilogy(rslt[B4]['upperstateenergy'], rslt[B4]['T_B'], 's', markerfacecolor='none', label='Band 4')
    ax1.semilogy(rslt[B7]['upperstateenergy'], rslt[B7]['T_B'], 's', markerfacecolor='none', label='Band 7')
    ax1.set_xlabel("E$_U$ [K]")
    ax1.set_ylabel("T$_B$ [K]")
    ax1.set_ylim(0.1, 200)
    pl.legend(loc='best')
    pl.tight_layout()

    ax2 = pl.subplot(2,1,2)
    L0, = ax2.semilogy(rslt[vzero]['upperstateenergy'], rslt[vzero]['upperlevelpop'], ',')
    L1, = ax2.semilogy(rslt[vone]['upperstateenergy'], rslt[vone]['upperlevelpop'], ',')
    L2, =ax2.semilogy(rslt[vtwo]['upperstateenergy'], rslt[vtwo]['upperlevelpop'], ',')
    ax2.semilogy(rslt[vzero & obs]['upperstateenergy'], rslt[vzero & obs]['upperlevelpop'], 'o', color=L0.get_color(), label='v=0')
    ax2.semilogy(rslt[vone & obs]['upperstateenergy'], rslt[vone & obs]['upperlevelpop'], 'o', color=L1.get_color(), label='v=1')
    ax2.semilogy(rslt[vtwo & obs]['upperstateenergy'], rslt[vtwo & obs]['upperlevelpop'], 'o', color=L2.get_color(), label='v=2')
    ax2.semilogy(rslt[B4]['upperstateenergy'], rslt[B4]['upperlevelpop'], 's', markerfacecolor='none', label='Band 4')
    ax2.semilogy(rslt[B7]['upperstateenergy'], rslt[B7]['upperlevelpop'], 's', markerfacecolor='none', label='Band 7')
    #ax2.semilogy(rslt[Jeight]['upperstateenergy'], rslt[Jeight]['upperlevelpop'], 's', markerfacecolor='none', label='J=8')
    ax2.set_xlabel("E$_U$ [K]")
    ax2.set_ylabel("N$_U$")
    #pl.ylim(0.1, 200)
    #pl.legend(loc='best')
    pl.tight_layout()

rr.background_brightness = artificial_29um
rslt = (rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=150*u.K, tbg=None))

pl.figure()
pl.title("Model 1: 29$\\mu$m bump")
makeplot(rslt)

rslt = (rr(density=1e10*u.cm**-3, column=2e14*u.cm**-2, temperature=150*u.K, tbg=2.73))

pl.figure()
pl.title("Model 2: High density, low-column")
makeplot(rslt)

rslt = (rr(density=1e10*u.cm**-3, column=1e18*u.cm**-2, temperature=150*u.K, tbg=2.73))

pl.figure()
pl.title("Model 3: Optically Thick, high-density high-column")
makeplot(rslt)


rslt = (rr(density=1e11*u.cm**-3, column=1e19*u.cm**-2, temperature=150*u.K, tbg=2.73))

pl.figure()
pl.title("Model 4: Optically Thicker, higher-density higher-column")
makeplot(rslt)

rslt = (rr(density=1e4*u.cm**-3, column=1e20*u.cm**-2, temperature=150*u.K, tbg=2.73))

pl.figure()
pl.title("Model 5: Low-density, High-column")
makeplot(rslt)

artificial_295um = bgfunc(freq)
artificial_295um[rovib_range] += 1e-9*rr.background_brightness.unit * np.exp(-(wl[rovib_range]-29.5*u.um)**2/(2*(2.5*u.um)**2))
artificial_285um = bgfunc(freq)
artificial_285um[rovib_range] += 1e-9*rr.background_brightness.unit * np.exp(-(wl[rovib_range]-28.5*u.um)**2/(2*(2.5*u.um)**2))


rr.background_brightness = artificial_285um
rslt = (rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=150*u.K, tbg=None))
pl.figure()
pl.title("Model 6: 28.5$\\mu$m bump")
makeplot(rslt)

rr.background_brightness = artificial_295um
rslt = (rr(density=1e4*u.cm**-3, column=2e14*u.cm**-2, temperature=150*u.K, tbg=None))
pl.figure()
pl.title("Model 7: 29.5$\\mu$m bump")
makeplot(rslt)

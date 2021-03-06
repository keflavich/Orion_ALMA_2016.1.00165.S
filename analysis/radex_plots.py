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

for density in (1e2, 1e5, 1e8, 1e11, 1e14):
    rslt = rr(density={'H2':density})
    critdens = {rr.quantum_number[iupp]:
                rr.radex.rmolec.aeinst[rr.radex.imolec.iupp == iupp].sum() / (rr.radex.collie.ctot[iupp-1] / rr.total_density)
                for iupp in np.arange(1, rr.radex.imolec.iupp.max())}
    print(density, critdens[b'0_7   '], rr.level_population[7], critdens[b'1_7   '], rr.level_population[68], "\n", rslt[v0_76 | v1_76])


v0_1817 = (rslt['upperlevel'] == '0_18  ') & (rslt['lowerlevel'] == '0_17  ')
v1_1817 = (rslt['upperlevel'] == '1_18  ') & (rslt['lowerlevel'] == '1_17  ')
v2_1817 = (rslt['upperlevel'] == '2_18  ') & (rslt['lowerlevel'] == '2_17  ')
v3_1817 = (rslt['upperlevel'] == '3_18  ') & (rslt['lowerlevel'] == '3_17  ')

obs = (((85.5 < rslt['frequency']) & (89.5 > rslt['frequency'])) |
       ((97.5 < rslt['frequency']) & (101.5 > rslt['frequency'])) |
       ((229.0 < rslt['frequency']) & (233.7 > rslt['frequency'])) |
       ((214.25 < rslt['frequency']) & (218.9 > rslt['frequency'])) |
       ((344.1 < rslt['frequency']) & (348 > rslt['frequency'])) |
       ((332.1 < rslt['frequency']) & (335.8 > rslt['frequency'])))


# play with different backgrounds

print(rr(density=1e4*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=2.73)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
print(chi2(rr.get_table()))

freq = rr.frequency
wl = rr.frequency.to(u.um, u.spectral())
bb100_half_plus_cmb = dust_emissivity.blackbody.blackbody(nu=freq, temperature=100*u.K)/2. + dust_emissivity.blackbody.blackbody(nu=freq, temperature=2.73*u.K)

rr.background_brightness = bb100_half_plus_cmb
print(rr(density=1e4*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=None)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
print(chi2(rr.get_table()))

bb4000smallff_100_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=100*u.K)/2. +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=2.73*u.K))

rr.background_brightness = bb4000smallff_100_half_plus_cmb
print(rr(density=1e4*u.cm**-3, column=1e14*u.cm**-2, temperature=100*u.K, tbg=None)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
print(chi2(rr.get_table()))

bb4000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=200*u.K)/2. +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                   dust_emissivity.blackbody.blackbody(nu=freq,
                                                                       temperature=2.73*u.K))
rr.background_brightness = bb4000smallff_200_half_plus_cmb
print(rr(density=1e5*u.cm**-3, column=1e15*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))
print(rr(density=1e6*u.cm**-3, column=1e15*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))
print(rr(density=1e7*u.cm**-3, column=1e15*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))
print(rr(density=1e5*u.cm**-3, column=1e16*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))
print(rr(density=1e6*u.cm**-3, column=1e16*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))
print(rr(density=1e7*u.cm**-3, column=1e16*u.cm**-2, temperature=100*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))

mbb1000smallff_200_half_plus_cmb = (dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=200*u.K)/2. +
                                    dust_emissivity.blackbody.modified_blackbody(nu=freq,
                                                                                 temperature=1000*u.K)*(10*u.au)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=4000*u.K)*(100*u.R_sun)**2/(30*u.au)**2 +
                                    dust_emissivity.blackbody.blackbody(nu=freq,
                                                                        temperature=2.73*u.K))
rr.background_brightness = mbb1000smallff_200_half_plus_cmb
print(rr(density=5e5*u.cm**-3, column=1e15*u.cm**-2, temperature=150*u.K, tbg=None)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32])
print(chi2(rr.get_table()))
print(rr(density=5e5*u.cm**-3, column=1e15*u.cm**-2, temperature=150*u.K, tbg=None)[obs])
print(chi2(rr.get_table()))

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


if False:

    densities = np.logspace(5,13,50)

    if 'data_N1e12_T100' not in locals():
        data_N1e12_T100 = [rr(density=density,
                              column=1e12,
                              temperature=100)[
                                  v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32]
                           for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e12_T100]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e12_T100]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e12_T100]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e12_T100]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e12_T100]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e12_T100]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e12_T100]), label='v=3-2 J=7-6')
    pl.ylim(0,110)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e12_T100]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e12_T100]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e12_T100]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e12_T100]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_N=1e12.png'))


    densities = np.logspace(5,13,50)

    if 'data_N1e13_T100' not in locals():
        data_N1e13_T100 = [rr(density=density, column=1e13, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

    pl.clf()
    pl.subplot(3,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e13_T100]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e13_T100]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e13_T100]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e13_T100]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e13_T100]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e13_T100]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e13_T100]), label='v=3-2 J=7-6')
    pl.ylim(0,110)
    pl.legend(loc='lower right')

    pl.subplot(3,1,2)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{B}$ [K]")
    pl.semilogx(densities, np.array([x[3]['T_B'] for x in data_N1e13_T100]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['T_B'] for x in data_N1e13_T100]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['T_B'] for x in data_N1e13_T100]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['T_B'] for x in data_N1e13_T100]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['T_B'] for x in data_N1e13_T100]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['T_B'] for x in data_N1e13_T100]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['T_B'] for x in data_N1e13_T100]), label='v=3-2 J=7-6')
    pl.ylim(0,110)
    pl.legend(loc='lower right')


    pl.subplot(3,1,3)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e13_T100]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e13_T100]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e13_T100]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e13_T100]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_N=1e13.png'))



    if 'data_N1e12_T1000' not in locals():
        data_N1e12_T1000 = [rr(density=density, column=1e12, temperature=1000)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e12_T1000]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e12_T1000]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e12_T1000]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e12_T1000]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e12_T1000]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e12_T1000]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e12_T1000]), label='v=3-2 J=7-6')
    pl.ylim(0,1100)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e12_T1000]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e12_T1000]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e12_T1000]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e12_T1000]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_1000K_N=1e12.png'))



    densities = np.logspace(5,13,50)

    if 'data_N1e15_T100' not in locals():
        data_N1e15_T100 = [rr(density=density, column=1e15, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e15_T100]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e15_T100]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e15_T100]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e15_T100]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e15_T100]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e15_T100]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e15_T100]), label='v=3-2 J=7-6')
    pl.ylim(0,110)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e15_T100]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e15_T100]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e15_T100]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e15_T100]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_N=1e15.png'))


    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    densities = np.logspace(5,13,50)

    if 'data_N1e14_T100' not in locals():
        data_N1e14_T100 = [rr(density=density, column=1e14, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e14_T100]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e14_T100]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e14_T100]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e14_T100]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e14_T100]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e12_T100]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e12_T100]), label='v=3-2 J=7-6')
    pl.ylim(0,110)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e14_T100]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e14_T100]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e14_T100]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e14_T100]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_N=1e14.png'))



    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    if 'data_N1e12_T100_J1817' not in locals():
        data_N1e12_T100_J1817 = [rr(density=density, column=1e12, temperature=100)[v0_1817 | v1_1817 | v2_1817 | v3_1817] for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_B$ [K]")
    pl.semilogx(densities, np.array([x[3]['T_B'] for x in data_N1e12_T100_J1817]), label='v=0 J=18-17')
    pl.semilogx(densities, np.array([x[2]['T_B'] for x in data_N1e12_T100_J1817]), label='v=1 J=18-17')
    pl.semilogx(densities, np.array([x[1]['T_B'] for x in data_N1e12_T100_J1817]), label='v=2 J=18-17')
    pl.semilogx(densities, np.array([x[0]['T_B'] for x in data_N1e12_T100_J1817]), label='v=3 J=18-17')
    pl.ylim(0,10)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e12_T100_J1817]), label='v=0 J=18-17')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e12_T100_J1817]), label='v=1 J=18-17')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e12_T100_J1817]), label='v=2 J=18-17')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e12_T100_J1817]), label='v=3 J=18-17')
    pl.ylim(0,110)
    pl.legend(loc='lower right')
    pl.savefig(paths.fpath('radex/NaCl_J=18-17_vs_density_100K_N=1e12.png'))




    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    if 'data_N1e13_T100_J1817' not in locals():
        data_N1e13_T100_J1817 = [rr(density=density, column=1e13, temperature=100)[v0_1817 | v1_1817 | v2_1817 | v3_1817] for density in densities]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_B$ [K]")
    pl.semilogx(densities, np.array([x[3]['T_B'] for x in data_N1e13_T100_J1817]), label='v=0 J=18-17')
    pl.semilogx(densities, np.array([x[2]['T_B'] for x in data_N1e13_T100_J1817]), label='v=1 J=18-17')
    pl.semilogx(densities, np.array([x[1]['T_B'] for x in data_N1e13_T100_J1817]), label='v=2 J=18-17')
    pl.semilogx(densities, np.array([x[0]['T_B'] for x in data_N1e13_T100_J1817]), label='v=3 J=18-17')
    pl.ylim(0,40)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e13_T100_J1817]), label='v=0 J=18-17')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e13_T100_J1817]), label='v=1 J=18-17')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e13_T100_J1817]), label='v=2 J=18-17')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e13_T100_J1817]), label='v=3 J=18-17')
    pl.ylim(0,110)
    pl.legend(loc='lower right')
    pl.savefig(paths.fpath('radex/NaCl_J=18-17_vs_density_100K_N=1e13.png'))


    columns = np.logspace(12, 20, 50)

    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    if 'data_n1e8_T100' not in locals():
        data_n1e8_T100 = [rr(density=1e8, column=column, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for column in columns]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Column of NaCl [cm$^{-2}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(columns, np.array([x[3]['Tex'] for x in data_n1e8_T100]), label='v=0 J=7-6')
    pl.semilogx(columns, np.array([x[2]['Tex'] for x in data_n1e8_T100]), label='v=1 J=7-6')
    pl.semilogx(columns, np.array([x[1]['Tex'] for x in data_n1e8_T100]), label='v=2 J=7-6')
    pl.semilogx(columns, np.array([x[0]['Tex'] for x in data_n1e8_T100]), label='v=3 J=7-6')
    pl.semilogx(columns, np.array([x[6]['Tex'] for x in data_n1e8_T100]), label='v=1-0 J=7-6')
    pl.semilogx(columns, np.array([x[5]['Tex'] for x in data_n1e8_T100]), label='v=2-1 J=7-6')
    pl.semilogx(columns, np.array([x[4]['Tex'] for x in data_n1e8_T100]), label='v=3-2 J=7-6')
    pl.ylim(-5,110)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, columns.min(), columns.max(), color='k', linestyle='--')
    pl.loglog(columns, np.array([x[3]['tau'] for x in data_n1e8_T100]), label='v=0 J=7-6')
    pl.loglog(columns, np.array([x[4]['tau'] for x in data_n1e8_T100]), label='v=3-2 J=7-6')
    pl.loglog(columns, np.array([x[5]['tau'] for x in data_n1e8_T100]), label='v=2-1 J=7-6')
    pl.loglog(columns, np.array([x[6]['tau'] for x in data_n1e8_T100]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Column of NaCl [cm$^{-2}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_column_100K_n=1e8.png'))

    columns = np.logspace(12, 20, 50)

    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    if 'data_n1e8_T1000' not in locals():
        data_n1e8_T1000 = [rr(density=1e8, column=column, temperature=1000)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for column in columns]

    pl.clf()
    pl.subplot(2,1,1)
    pl.xlabel("Column of NaCl [cm$^{-2}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(columns, np.array([x[3]['Tex'] for x in data_n1e8_T1000]), label='v=0 J=7-6')
    pl.semilogx(columns, np.array([x[2]['Tex'] for x in data_n1e8_T1000]), label='v=1 J=7-6')
    pl.semilogx(columns, np.array([x[1]['Tex'] for x in data_n1e8_T1000]), label='v=2 J=7-6')
    pl.semilogx(columns, np.array([x[0]['Tex'] for x in data_n1e8_T1000]), label='v=3 J=7-6')
    pl.semilogx(columns, np.array([x[6]['Tex'] for x in data_n1e8_T1000]), label='v=1-0 J=7-6')
    pl.semilogx(columns, np.array([x[5]['Tex'] for x in data_n1e8_T1000]), label='v=2-1 J=7-6')
    pl.semilogx(columns, np.array([x[4]['Tex'] for x in data_n1e8_T1000]), label='v=3-2 J=7-6')
    pl.ylim(-5,110)
    pl.legend(loc='lower right')
    pl.subplot(2,1,2)
    pl.hlines(1, columns.min(), columns.max(), color='k', linestyle='--')
    pl.loglog(columns, np.array([x[3]['tau'] for x in data_n1e8_T1000]), label='v=0 J=7-6')
    pl.loglog(columns, np.array([x[4]['tau'] for x in data_n1e8_T1000]), label='v=3-2 J=7-6')
    pl.loglog(columns, np.array([x[5]['tau'] for x in data_n1e8_T1000]), label='v=2-1 J=7-6')
    pl.loglog(columns, np.array([x[6]['tau'] for x in data_n1e8_T1000]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right')
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Column of NaCl [cm$^{-2}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_column_1000K_n=1e8.png'))

    densities = np.logspace(3,13,50)

    #rr = lambda: pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=1e14)
    if 'data_N1e14_T1000' not in locals():
        data_N1e14_T1000 = [rr(density={'H2':density}, column=1e14, temperature=1000)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

    pl.clf()
    pl.subplot(3,1,1)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{ex}$ [K]")
    pl.semilogx(densities, np.array([x[3]['Tex'] for x in data_N1e14_T1000]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['Tex'] for x in data_N1e14_T1000]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['Tex'] for x in data_N1e14_T1000]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['Tex'] for x in data_N1e14_T1000]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['Tex'] for x in data_N1e14_T1000]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['Tex'] for x in data_N1e14_T1000]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['Tex'] for x in data_N1e14_T1000]), label='v=3-2 J=7-6')
    pl.ylim(0,1100)
    pl.legend(loc='lower right', fontsize=8)

    pl.subplot(3,1,2)
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.ylabel("T$_{B}$ [K]")
    pl.semilogx(densities, np.array([x[3]['T_B'] for x in data_N1e14_T1000]), label='v=0 J=7-6')
    pl.semilogx(densities, np.array([x[2]['T_B'] for x in data_N1e14_T1000]), label='v=1 J=7-6')
    pl.semilogx(densities, np.array([x[1]['T_B'] for x in data_N1e14_T1000]), label='v=2 J=7-6')
    pl.semilogx(densities, np.array([x[0]['T_B'] for x in data_N1e14_T1000]), label='v=3 J=7-6')
    pl.semilogx(densities, np.array([x[6]['T_B'] for x in data_N1e14_T1000]), label='v=1-0 J=7-6')
    pl.semilogx(densities, np.array([x[5]['T_B'] for x in data_N1e14_T1000]), label='v=2-1 J=7-6')
    pl.semilogx(densities, np.array([x[4]['T_B'] for x in data_N1e14_T1000]), label='v=3-2 J=7-6')
    pl.ylim(0,10)
    pl.legend(loc='lower right', fontsize=8)


    pl.subplot(3,1,3)
    pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
    pl.loglog(densities, np.array([x[3]['tau'] for x in data_N1e14_T1000]), label='v=0 J=7-6')
    pl.loglog(densities, np.array([x[4]['tau'] for x in data_N1e14_T1000]), label='v=3-2 J=7-6')
    pl.loglog(densities, np.array([x[5]['tau'] for x in data_N1e14_T1000]), label='v=2-1 J=7-6')
    pl.loglog(densities, np.array([x[6]['tau'] for x in data_N1e14_T1000]), label='v=1-0 J=7-6')
    pl.legend(loc='lower right', fontsize=8)
    pl.ylim(1e-7, 1e2)
    pl.ylabel("Optical Depth")
    pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
    pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_1000K_N=1e14.png'))



    # different sorts of plots....


    for ii,temperature in enumerate((200, 500, 1000)):
        pl.figure(ii)
        pl.clf()
        tbl = Table.read(paths.tpath('line_fits.txt'), format='ascii.fixed_width')
        naclmask = np.array(['23Na-35Cl' in x for x in tbl['Line Name']])
        pl.plot(tbl['Frequency'][naclmask], tbl['Fitted Amplitude K'][naclmask], 's', markerfacecolor='none', markeredgecolor='k')
        for (density, temperature, column) in (
            (1e4, temperature, 1e14),
            (1e9, temperature, 1e14),
            (1e10, temperature, 1e14),
            (1e11, temperature, 1e14),
            (1e12, temperature, 1e14),
            ):

            rslt = rr(density={'H2':density}, temperature=temperature, column=column)
            pl.plot(rslt['frequency'], rslt['T_B'], '.', label="n={0} T={1} N={2}".format(np.log10(density), temperature, np.log10(column)))
        pl.legend(loc='best')
        pl.xlim(0, 600)
        pl.ylim(-10, 200)
        pl.xlabel("Frequency")
        pl.ylabel("T$_B$")


    rr = pyradex.Radex(species='nacl', temperature=1000, density=1e8, column=4e13)
    bins = np.linspace(-18, -9.5)
    pl.figure(4)
    pl.clf()
    pl.title("T=1000 K")
    pl.hist(np.log10(rr.radex.collie.crate[:61,:61]/rr.total_density.value).ravel(), bins=bins, label='v=0')
    pl.hist(np.log10(rr.radex.collie.crate[:61,61:122]/rr.total_density.value).ravel(), bins=bins, alpha=0.25, label='v=0->1')
    pl.hist(np.log10(rr.radex.collie.crate[61:122,:61]/rr.total_density.value).ravel(), bins=bins, alpha=0.25, label='v=1->0')
    pl.hist(np.log10(rr.radex.collie.crate[:61,122:183]/rr.total_density.value).ravel(), bins=np.linspace(-15,-8), alpha=0.25, label='v=0->2')
    pl.hist(np.log10(rr.radex.collie.crate[122:183,:61]/rr.total_density.value).ravel(), bins=np.linspace(-15,-8), alpha=0.25, label='v=2->0')
    pl.legend(loc='best')
    pl.xlabel("log(A$_{ij}$")

    rr = pyradex.Radex(species='nacl', temperature=100, density=1e8, column=4e13)
    pl.figure(5)
    pl.clf()
    pl.title("T=100 K")
    pl.hist(np.log10(rr.radex.collie.crate[:61,:61]/rr.total_density.value).ravel(), bins=bins, label='v=0')
    pl.hist(np.log10(rr.radex.collie.crate[:61,61:122]/rr.total_density.value).ravel(), bins=bins, alpha=0.25, label='v=0->1')
    pl.hist(np.log10(rr.radex.collie.crate[61:122,:61]/rr.total_density.value).ravel(), bins=bins, alpha=0.25, label='v=1->0')
    pl.hist(np.log10(rr.radex.collie.crate[:61,122:183]/rr.total_density.value).ravel(), bins=np.linspace(-15,-8), alpha=0.25, label='v=0->2')
    pl.hist(np.log10(rr.radex.collie.crate[122:183,:61]/rr.total_density.value).ravel(), bins=np.linspace(-15,-8), alpha=0.25, label='v=2->0')
    pl.legend(loc='best')
    pl.xlabel("log(A$_{ij}$")

import pyradex
import pylab as pl
import numpy as np
import paths

rr = pyradex.Radex(species='nacl', temperature=100, density=1e8, abundance=1e-10)
rslt = rr()

v0_76 = (rslt['upperlevel'] == '0_7   ') & (rslt['lowerlevel'] == '0_6   ')
v1_76 = (rslt['upperlevel'] == '1_7   ') & (rslt['lowerlevel'] == '1_6   ')
v2_76 = (rslt['upperlevel'] == '2_7   ') & (rslt['lowerlevel'] == '2_6   ')
v3_76 = (rslt['upperlevel'] == '3_7   ') & (rslt['lowerlevel'] == '3_6   ')
v10 = (rslt['upperlevel'] == '1_1   ') & (rslt['lowerlevel'] == '0_0   ')
v21 = (rslt['upperlevel'] == '2_1   ') & (rslt['lowerlevel'] == '1_0   ')
v32 = (rslt['upperlevel'] == '3_1   ') & (rslt['lowerlevel'] == '2_0   ')

v0_1817 = (rslt['upperlevel'] == '0_18  ') & (rslt['lowerlevel'] == '0_17  ')
v1_1817 = (rslt['upperlevel'] == '1_18  ') & (rslt['lowerlevel'] == '1_17  ')
v2_1817 = (rslt['upperlevel'] == '2_18  ') & (rslt['lowerlevel'] == '2_17  ')
v3_1817 = (rslt['upperlevel'] == '3_18  ') & (rslt['lowerlevel'] == '3_17  ')


densities = np.logspace(6,13,100)

data = [rr(density=density, abundance=1e-10, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

pl.clf()
pl.subplot(2,1,1)
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=7-6')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=7-6')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=7-6')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=7-6')
pl.semilogx(densities, np.array([x[6]['T_B'] for x in data]), label='v=1-0')
pl.ylim(0,110)
pl.legend(loc='best')
pl.subplot(2,1,2)
pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
pl.loglog(densities, np.array([x[3]['tau'] for x in data]), label='v=0 J=7-6')
pl.loglog(densities, np.array([x[4]['tau'] for x in data]), label='v=3-2')
pl.loglog(densities, np.array([x[5]['tau'] for x in data]), label='v=2-1')
pl.loglog(densities, np.array([x[6]['tau'] for x in data]), label='v=1-0')
pl.legend(loc='best')
pl.ylim(1e-3, 1e5)
pl.ylabel("Optical Depth of $\Delta v = 1$")
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K.png'))


densities = np.logspace(6,13,100)

data = [rr(density=density, abundance=1e-9, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

pl.clf()
pl.subplot(2,1,1)
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=7-6')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=7-6')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=7-6')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=7-6')
pl.semilogx(densities, np.array([x[6]['T_B'] for x in data]), label='v=1-0')
pl.ylim(0,110)
pl.legend(loc='best')
pl.subplot(2,1,2)
pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
pl.loglog(densities, np.array([x[3]['tau'] for x in data]), label='v=0 J=7-6')
pl.loglog(densities, np.array([x[4]['tau'] for x in data]), label='v=3-2')
pl.loglog(densities, np.array([x[5]['tau'] for x in data]), label='v=2-1')
pl.loglog(densities, np.array([x[6]['tau'] for x in data]), label='v=1-0')
pl.legend(loc='best')
pl.ylim(1e-3, 1e5)
pl.ylabel("Optical Depth of $\Delta v = 1$")
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_x-9.png'))



data = [rr(density=density, abundance=1e-10, temperature=1000)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

pl.clf()
pl.subplot(2,1,1)
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=7-6')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=7-6')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=7-6')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=7-6')
pl.semilogx(densities, np.array([x[6]['T_B'] for x in data]), label='v=1-0')
pl.ylim(0,1100)
pl.legend(loc='best')
pl.subplot(2,1,2)
pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
pl.loglog(densities, np.array([x[3]['tau'] for x in data]), label='v=0 J=7-6')
pl.loglog(densities, np.array([x[4]['tau'] for x in data]), label='v=3-2')
pl.loglog(densities, np.array([x[5]['tau'] for x in data]), label='v=2-1')
pl.loglog(densities, np.array([x[6]['tau'] for x in data]), label='v=1-0')
pl.legend(loc='best')
pl.ylim(1e-3, 1e5)
pl.ylabel("Optical Depth of $\Delta v = 1$")
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_1000K.png'))



densities = np.logspace(6,13,100)

data = [rr(density=density, abundance=1e-7, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

pl.clf()
pl.subplot(2,1,1)
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=7-6')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=7-6')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=7-6')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=7-6')
pl.semilogx(densities, np.array([x[6]['T_B'] for x in data]), label='v=1-0')
pl.ylim(0,110)
pl.legend(loc='best')
pl.subplot(2,1,2)
pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
pl.loglog(densities, np.array([x[3]['tau'] for x in data]), label='v=0 J=7-6')
pl.loglog(densities, np.array([x[4]['tau'] for x in data]), label='v=3-2')
pl.loglog(densities, np.array([x[5]['tau'] for x in data]), label='v=2-1')
pl.loglog(densities, np.array([x[6]['tau'] for x in data]), label='v=1-0')
pl.legend(loc='best')
pl.ylim(1e-3, 1e5)
pl.ylabel("Optical Depth of $\Delta v = 1$")
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_X-7.png'))


densities = np.logspace(6,13,100)

data = [rr(density=density, abundance=1e-8, temperature=100)[v0_76 | v1_76 | v2_76 | v3_76 | v10 | v21 | v32] for density in densities]

pl.clf()
pl.subplot(2,1,1)
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=7-6')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=7-6')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=7-6')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=7-6')
pl.semilogx(densities, np.array([x[6]['T_B'] for x in data]), label='v=1-0')
pl.ylim(0,110)
pl.legend(loc='best')
pl.subplot(2,1,2)
pl.hlines(1, 1e6, 1e13, color='k', linestyle='--')
pl.loglog(densities, np.array([x[3]['tau'] for x in data]), label='v=0 J=7-6')
pl.loglog(densities, np.array([x[4]['tau'] for x in data]), label='v=3-2')
pl.loglog(densities, np.array([x[5]['tau'] for x in data]), label='v=2-1')
pl.loglog(densities, np.array([x[6]['tau'] for x in data]), label='v=1-0')
pl.legend(loc='best')
pl.ylim(1e-3, 1e5)
pl.ylabel("Optical Depth of $\Delta v = 1$")
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.savefig(paths.fpath('radex/NaCl_J=7-6_vs_density_100K_X-8.png'))



data = [rr(density=density, abundance=1e-10, temperature=100)[v0_1817 | v1_1817 | v2_1817 | v3_1817] for density in densities]

pl.clf()
pl.xlabel("Density of H$_2$ [cm$^{-3}$]")
pl.ylabel("T$_B$ [K]")
pl.semilogx(densities, np.array([x[3]['T_B'] for x in data]), label='v=0 J=18-17')
pl.semilogx(densities, np.array([x[2]['T_B'] for x in data]), label='v=1 J=18-17')
pl.semilogx(densities, np.array([x[1]['T_B'] for x in data]), label='v=2 J=18-17')
pl.semilogx(densities, np.array([x[0]['T_B'] for x in data]), label='v=3 J=18-17')
pl.ylim(0,110)
pl.legend(loc='best')
pl.savefig(paths.fpath('radex/NaCl_J=18-17_vs_density_100K.png'))



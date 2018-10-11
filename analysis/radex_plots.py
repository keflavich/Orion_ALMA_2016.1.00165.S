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
v10 = (rslt['upperlevel'] == '1_7   ') & (rslt['lowerlevel'] == '0_6   ')
v21 = (rslt['upperlevel'] == '2_7   ') & (rslt['lowerlevel'] == '1_6   ')
v32 = (rslt['upperlevel'] == '3_7   ') & (rslt['lowerlevel'] == '2_6   ')

v0_1817 = (rslt['upperlevel'] == '0_18  ') & (rslt['lowerlevel'] == '0_17  ')
v1_1817 = (rslt['upperlevel'] == '1_18  ') & (rslt['lowerlevel'] == '1_17  ')
v2_1817 = (rslt['upperlevel'] == '2_18  ') & (rslt['lowerlevel'] == '2_17  ')
v3_1817 = (rslt['upperlevel'] == '3_18  ') & (rslt['lowerlevel'] == '3_17  ')


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
pl.subplot(2,1,1)
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
pl.subplot(2,1,2)
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

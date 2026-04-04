import numpy as np
from astropy import constants
from astropy import table
from astropy import units as u

import paths
from lte_modeling_tools import (rovib_lte_model_generator,
                                simple_lte_model_generator, kkms_of_nupper,
                                nupper_of_kkms)


tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

kcl35mask = np.array([(not hasattr(row['Species'], 'mask')) and
                     ('KCl' == row['Species'][:3] or
                      '39K-35Cl' in row['Species']) for row in tbl])
kcl35tbl = tbl[kcl35mask]

v0 = np.array(['v=0' in row['Species'] for row in kcl35tbl])


# Using LTE modeling tools, what integrated intensity do we expect for each
# line for some arbitrary column density (we are ignoring optical depth) for a
# high excitation temperature?
#
# this is split out from kcl_rotational_diagrams.py
mod = simple_lte_model_generator()
mod.tem = 1000
mod.logcolumn = np.log(1e9)
fluxes = kkms_of_nupper(np.exp(mod(kcl35tbl[v0]['EU_K'])),
                        kcl35tbl[v0]['Frequency'],
                        kcl35tbl[v0]['Aij'],
                        #kcl35tbl[v0]['deg'],
                        )
#Out[363]: <Quantity [ 12.87528571,  53.19600633,  59.7568103 , 109.06635309] K km / s>
s345 = fluxes[-1]
s100 = fluxes[0]
print("For T_ex = 1000 K, s100 / s345 = {0}".format(s100/s345))
print("                   s345 / s100 = {0}".format(s345/s100))

# the flux measurements at 100 GHz in the KCl 13-12 and 345 GHz 45-44 lines, in
# Kelvin units
f100 = 15.3
f345 = 9.7

for beta in (1., 1.5, 2):
    print()
    print()
    print(f"          beta={beta}")

    # this is the opacity ratio between the two frequencies for a given dust Beta
    # kappa = kappa_0 * nu^beta
    ratio = (345/100)**beta

    # eqn is: F_1 / F_2 = S_1 exp(-kappa1 N) / (S_2 exp(-kappa2 N)
    # F1/F2 * S_2/S_1 = exp(-kappa1 N + kappa2 N)
    # kappa1 = kappa2 * ratio
    # log F1S2/F2S1 = N (kappa2 * (1-ratio))
    # log F1S2/F2S1 = N (kappa1 * (1/ratio-1))
    # assume S_1 = S_2

    logf1f2 = np.log(f345/f100)
    tau345 = kappa1N = logf1f2 / (1/ratio - 1)
    tau100 = kappa2N = logf1f2 / (1 - ratio)

    print()
    print("Optical depths for F100 = F345")
    print(f"tau345 = {tau345}")
    print(f"tau100 = {tau100}")
    print("tau345/tau100 = {0}".format(tau345/tau100))



    # these are th evalues pulled from the T=1000K analysis above
    s345 = 109.06
    s100 = 12.875

    logf1f2 = np.log(f345/f100*s100/s345)
    tau345 = kappa1N = logf1f2 / (1/ratio - 1)
    tau100 = kappa2N = logf1f2 / (1 - ratio)

    print()
    print("Optical depths for F100 = {0}F345, i.e., for LTE T=1000 K".format(s100/s345))
    print(f"tau345 = {tau345}")
    print(f"tau100 = {tau100}")
    print("tau345/tau100 = {0}".format(tau345/tau100))


    print("\n\n2025-10-17: What if we don't ignore line optical depth?")
    print("For Orion-like cases, the lowest 3 vibrational states are optically thick - which explains their constant(ish) observed brightness")

    from pyspeckit.spectrum.models.lte_molecule import line_tau, get_molecular_parameters

    import pylab as pl
    pl.figure(figsize=(14,6))
    pl.clf()

    dV = 20 * u.km/u.s
    for dV in (7, 20)*u.km/u.s:
        total_column = 1e17*u.cm**-2
        pl.clf()
        for ii, tex in enumerate([100, 250, 1000]*u.K):
            #freqs, aij, deg, EU, partfunc, jpltbl = get_molecular_parameters('NaCl', tex=tex)
            frqs, einsteinAij, degeneracies, EU, partfunc, cdmstbl = get_molecular_parameters('NaCl, v=0-15', catalog='CDMS', fmin=5*u.GHz, fmax=360*u.GHz, return_table=True)

            # we want to go from integral(tau dnu) -> tau(max)
            # using 20 km/s to match Kei's values
            dnu = (dV / constants.c) * frqs

            EU = u.Quantity(EU, u.erg)
            einsteinAij = u.Quantity(10**einsteinAij, 1/u.s)
            taus = line_tau(tex=tex, total_column=total_column, partition_function=partfunc(tex), degeneracy=degeneracies, frequency=frqs,
                            energy_upper=EU, einstein_A=einsteinAij)
            cdmstbl['taus'] = taus
            pl.subplot(1, 3, ii+1)
            for vu in np.arange(9):
                sel = (cdmstbl['Ku'] == vu) & (cdmstbl['Kl'] == vu)
                pl.semilogy(frqs[sel].to(u.GHz), taus[sel] / dnu[sel].to(u.Hz).value, 'o', label=f'$v_u=${vu}');
            pl.axhline(1, color='k', linestyle='--')
            pl.title(f"T={tex.to_string(format='latex')}")
            pl.xlabel("Frequency")
            pl.ylabel(r"$\tau$")
        pl.suptitle(f"$N=${total_column.to_string(format='latex')}  and  dV = {dV.to_string(format='latex')}", )
        pl.legend(loc='upper left', bbox_to_anchor=(1,1));
        pl.tight_layout()
        pl.savefig(paths.fpath(f"opticaldepth_vs_frequency_N={np.log10(total_column.value)}_dV={dV.value}.png"), bbox_inches='tight')
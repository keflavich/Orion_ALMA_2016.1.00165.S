import numpy as np
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
                        kcl35tbl[v0]['deg'])
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

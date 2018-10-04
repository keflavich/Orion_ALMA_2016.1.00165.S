""" had to add this to the fortran file:
Cf2py intent(out) rr
to make f2py return the value 'rr'

then run
f2py -c -m rates_ios rates_ios.f
"""
import itertools
import numpy as np
from astropy.utils.console import ProgressBar
from astropy import units as u
from salt_tables import NaCl
import rates_ios

# can't have identical level energies
NaCl['E_U'][375] = 536.9248881998526
NaCl['E_U'][np.array([61,3670,3881])] = 536.5639967986373
NaCl['E_U'][np.array([278, 3304, 4292, 7166, 7771])] = 1339.31

# v=x-x J=1-0 energy level - frequency-to-kelvin
# (these are K values)
groundstates = {(0,0): 0,
                (1,0): 519.67949313,
                (2,0): 1034.28410924,
                (3,0): 1543.88869568,
               }

maxv = 3
maxj = 60
assert maxv*maxj < 99999 # max # collision rates

with open('nacl.dat', 'w') as fh:

    fh.write("""!MOLECULE
NaCl
!MOLECULAR WEIGHT
58
! NUMBER OF ENERGY LEVELS
{nlev}
!LEVEL + ENERGIES(cm^-1) + WEIGHT + J + V
""".format(nlev=(maxv + 1) * (maxj + 1)))


    ii = 1
    leveldict = {}
    levelenergy = {}
    for vv in range(0,maxv+1):
        if ii in leveldict.values():
            raise ValueError((vv,))
        leveldict[(vv,0)] = ii #1 + maxj*vv
        energy_K = groundstates[(vv,0)]
        energy = (energy_K*u.K).to(u.eV, u.temperature_energy()).to(u.cm**-1, u.spectral()).value
        levelenergy[(vv,0)] = energy
        degen = 16
        # manually write out the ground state, since it's not in the table
        fh.write("{0:5d} {1:15.9f} {2:5.1f} {3:4d}\n".format(ii, energy, degen, 0, vv))
        ii += 1
        for jj in range(1,maxj+1):
            row = NaCl[(NaCl['vu'] == vv) & (NaCl['Ju'] == jj)]
            energy = row['E_U'].quantity[0].to(u.eV, u.temperature_energy()).to(u.cm**-1, u.spectral()).value
            degen = 16 + 32 * jj
            fh.write("{0:5d} {1:15.9f} {2:5.1f} {3:4d}\n".format(ii, energy, degen, jj, vv))
            if ii in leveldict.values():
                raise ValueError((vv,jj))
            leveldict[(vv,jj)] = ii
            if energy in levelenergy.values():
                raise ValueError((vv,jj))
            levelenergy[(vv,jj)] = energy
            ii += 1

    radtrans_strings = []
    ii = 1
    for row in NaCl:
        ju,jl = row['Ju'], row['Jl']
        vu,vl = row['vu'], row['vl']
        if (vu, ju) in leveldict and (vl, jl) in leveldict:

            if levelenergy[(vu,ju)] - levelenergy[(vl,jl)] < 1e-4:
                raise ValueError

            # only include the levels we're interested in
            radtrans_strings.append(
                "{0:5d} {1:5d} {2:5d} {3:11.3e} {4:20.8f} {5:10.2f}\n"
                .format(ii, leveldict[(vu, ju)], leveldict[(vl, jl)],
                        row['Aij'], row['Freq'], row['E_U'])
            )
            ii += 1

    fh.write("""! NUMBER OF RADIATIVE TRANSITIONS
{0}
!TRANS + U + L + A(s^-1) + FREQ(GHz) + E_u/k(K)
""".format(len(radtrans_strings)))

    for rtstr in radtrans_strings:
        fh.write(rtstr)


    levels = [(v,j) for v in range(0,maxv+1) for j in range(0,maxj+1)]
    pairs = list(itertools.permutations(levels, 2))


    temperatures = (10,50,100,150,2000)

    def pos_or_zero(x):
        if x <= 0:
            return 1e-50
        else:
            return x

    pb = ProgressBar(len(pairs))

    ii = 1
    ratestrings = []
    for ((v1,j1), (v2,j2)) in (pairs):
        if (j1==j2):
            # delta-J = 0 are forbidden
            continue
        #sel1 = (NaCl['vu'] == v1) & (NaCl['Ju'] == j1)
        #sel2 = (NaCl['Ju'] == j2) & (NaCl['vu'] == v2)
        #if j2 == 0 and v2 == 0:
        #    energy2 = 0
        #else:
        #    energy2 = u.Quantity(NaCl[sel2]['E_U'][0], u.K).to(u.eV, u.temperature_energy())
        #energy1 = u.Quantity(NaCl[sel1]['E_U'], u.K).to(u.eV, u.temperature_energy())
        energy1 = levelenergy[(v1,j1)]
        energy2 = levelenergy[(v2,j2)]

        if energy1 - energy2 == 0:
            raise ValueError

        # why not both directions?
        # if v1<v2 and j1<=j2:
        #     # only include collisions in one direction
        #     # (v1 = vu)
        #     pb.update()
        #     continue
        #r100 = rates_ios.rates(v1, j1, v2, j2, 100)
        #rr100 = rates_ios.rates(v2, j2, v1, j1, 100)
        #if r100 < 1e-19 and rr100 < 1e-19:
        #    # Skip negligible rates to avoid singular matrices
        #    pb.update()
        #    continue
        rates = [pos_or_zero(rates_ios.rates(v1, j1, v2, j2, tem)) for tem in temperatures]
        ratestr = " ".join(["{0:7.1e}".format(rr) for rr in rates])
        ratestrings.append("{0:5d} {1:5d} {2:5d} {3}\n"
                           .format(ii, leveldict[(v1, j1)], leveldict[(v2, j2)], ratestr))
        ii += 1

        pb.update()

    fh.write("""!NUMBER OF COLL PARTNERS
1
!COLLISIONS BETWEEN
1 NaCl-H2 scaled from 1.38 * NaCl-He from rates_ios.f
!NUMBER OF COLL TRANS
{npairs}
!NUMBER OF COLL TEMPS
5
!COLL TEMPS
   10.0  50.0  100.0  150.0 2000.0
!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)
""".format(npairs=len(ratestrings)))

    for ratestr in ratestrings:
        fh.write(ratestr)

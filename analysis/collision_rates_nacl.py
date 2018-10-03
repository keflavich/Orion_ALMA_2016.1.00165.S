""" had to add this to the fortran file:
Cf2py intent(out) rr
to make f2py return the value 'rr'

then run
f2py -c -m rates_ios rates_ios.f
"""
import itertools
from astropy.utils.console import ProgressBar
from astropy import units as u
from salt_tables import NaCl
import rates_ios

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
""".format(nlev=(maxv+1) * maxj + 1))

    # manually write out the ground state, since it's not in the table
    fh.write("{0:5d} {1:15.9f} {2:5.1f} {3:4d}\n".format(1, 0, 16, 0, 0))

    ii = 2
    leveldict = {}
    for vv in range(0,maxv+1):
        for jj in range(1,maxj+1):
            row = NaCl[(NaCl['vu'] == vv) & (NaCl['Ju'] == jj)]
            energy = row['E_U'].quantity[0].to(u.eV, u.temperature_energy()).to(u.cm**-1, u.spectral()).value
            degen = 16 + 32 * jj
            fh.write("{0:5d} {1:15.9f} {2:5.1f} {3:4d}\n".format(ii, energy, degen, jj, vv))
            leveldict[(vv,jj)] = ii
            ii += 1

    radtrans_strings = []
    ii = 1
    for row in NaCl:
        ju,jl = row['Ju'], row['Jl']
        vu,vl = row['vu'], row['vl']
        if (vu, ju) in leveldict and (vl, jl) in leveldict:
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


    levels = [(v,j) for v in range(0,maxv+1) for j in range(1,maxj+1)]
    pairs = list(itertools.permutations(levels, 2))


    temperatures = list(map(float, "10.0   20.0   30.0   40.0   50.0   60.0   70.0   80.0   90.0  100.0  110.0  120.0  130.0  140.0  150.0  160.0  170.0  180.0  190.0  200.0  210.0  220.0  230.0  240.0  250.0  260.0  270.0  280.0  290.0  300.0  500.0 1000.0 2000.0".split()))

    def pos_or_zero(x):
        if x < 0:
            return 1e-60
        else:
            return x

    pb = ProgressBar(len(pairs))

    ii = 1
    ratestrings = []
    for ((v1,j1), (v2,j2)) in (pairs):
        r100 = rates_ios.rates(v1, j1, v2, j2, 100)
        if r100 < 1e-25:
            # Skip negligible rates to avoid singular matrices
            pb.update()
            continue
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
33
!COLL TEMPS
   10.0   20.0   30.0   40.0   50.0   60.0   70.0   80.0   90.0  100.0  110.0  120.0  130.0  140.0  150.0  160.0  170.0  180.0  190.0  200.0  210.0  220.0  230.0  240.0  250.0  260.0  270.0  280.0  290.0  300.0  500.0 1000.0 2000.0
!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)
""".format(npairs=len(ratestrings)))

    for ratestr in ratestrings:
        fh.write(ratestr)

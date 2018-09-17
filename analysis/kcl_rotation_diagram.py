import paths
import numpy as np
from astropy import units as u
from astropy import constants
from astropy import log
from astropy import table
from astropy import modeling

from astroquery.vamdc import Vamdc
from vamdclib import specmodel

import pylab as pl

import lines

def nupper_of_kkms(kkms, freq, Aul, degeneracies, replace_bad=None):
    """ Derived directly from pyspeckit eqns..."""

    if replace_bad:
        neg = kkms <= 0
        kkms[neg] = replace_bad

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)
    #nline = 1.95e3 * freq**2 / Aul * kkms
    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2
    # term2 = np.exp(-constants.h*freq/(constants.k_B*Tex)) -1
    # term2 -> kt / hnu
    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)
    return (nline * Khz / degeneracies).to(u.cm**-2)

def fit_tex(eupper, nupperoverg, verbose=False, plot=False, uplims=None,
            errors=None, min_nupper=1,
            replace_errors_with_uplims=False,
            molecule=None,
            color='r',
            marker='o',
            max_uplims='half',
            label='',
           ):
    """
    Fit the Boltzmann diagram

    Parameters
    ----------
    max_uplims: str or number
        The maximum number of upper limits before the fit is ignored completely
        and instead zeros are returned
    """
    model = modeling.models.Linear1D()
    #fitter = modeling.fitting.LevMarLSQFitter()
    fitter = modeling.fitting.LinearLSQFitter()

    nupperoverg_tofit = nupperoverg.copy().to(u.cm**-2).value

    if uplims is not None:
        upperlim_mask = nupperoverg < uplims

        # allow this magical keyword 'half'
        max_uplims = len(nupperoverg)/2. if max_uplims == 'half' else max_uplims

        if upperlim_mask.sum() > max_uplims:
            # too many upper limits = bad idea to fit.
            return 0*u.cm**-2, 0*u.K, 0, 0

        if errors is None:
            # if errors are not specified, we set the upper limits as values to
            # be fitted
            # (which gives a somewhat useful upper limit on the temperature)
            nupperoverg_tofit[upperlim_mask] = uplims[upperlim_mask]
        else:
            # otherwise, we set the values to zero-column and set the errors to
            # be whatever the upper limits are (hopefully a 1-sigma upper
            # limit)
            # 1.0 here becomes 0.0 in log and makes the relative errors meaningful
            nupperoverg_tofit[upperlim_mask] = 1.0
            if replace_errors_with_uplims:
                errors[upperlim_mask] = uplims[upperlim_mask]
    else:
        upperlim_mask = np.ones_like(nupperoverg_tofit, dtype='bool')

    # always ignore negatives & really low values
    good = nupperoverg_tofit > min_nupper
    # skip any fits that have fewer than 50% good values
    if good.sum() < len(nupperoverg_tofit)/2.:
        return 0*u.cm**-2, 0*u.K, 0, 0

    if errors is not None:
        rel_errors = errors / nupperoverg_tofit
        weights = 1. / rel_errors**2
        log.debug("Fitting with data = {0}, weights = {1}, errors = {2},"
                  "relative_errors = {3}"
                  .format(np.log(nupperoverg_tofit[good]),
                          np.log(weights[good]),
                          errors[good],
                          rel_errors[good],
                         ))
    else:
        # want log(weight) = 1
        weights = np.exp(np.ones_like(nupperoverg_tofit))

    result = fitter(model, eupper[good], np.log(nupperoverg_tofit[good]),
                    weights=np.log(weights[good]))
    tex = u.Quantity(-1./result.slope, u.K)

    partition_func = specmodel.calculate_partitionfunction(molecule.data['States'],
                                                           temperature=tex.value)
    assert len(partition_func) == 1
    Q_rot = tuple(partition_func.values())[0]

    Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2

    if verbose:
        print(("Tex={0}, Ntot={1}, log(Ntot)={4}, Q_rot={2}, "
               "nuplim={3}".format(tex, Ntot, Q_rot, upperlim_mask.sum(),
                                   np.log10(Ntot.value),
                                  )))

    if plot:
        import pylab as pl
        L, = pl.plot(eupper, np.log10(nupperoverg_tofit), marker=marker,
                     color=color, markeredgecolor='none', alpha=0.5,
                     linestyle='none',
                     #markersize=2,
                    )
        if uplims is not None:
            L, = pl.plot(eupper[upperlim_mask],
                         np.log10(uplims)[upperlim_mask], 'bv', alpha=0.2,
                         markersize=2)
            #L, = pl.plot(eupper[upperlim_mask],
            #             np.log10(nupperoverg)[upperlim_mask], 'bv', alpha=0.2)
        if errors is not None:
            yerr = np.array([np.log10(nupperoverg_tofit)-np.log10(nupperoverg_tofit-errors),
                             np.log10(nupperoverg_tofit+errors)-np.log10(nupperoverg_tofit)])
            # if lower limit is nan, set to zero
            yerr[0,:] = np.nan_to_num(yerr[0,:])
            if np.any(np.isnan(yerr[1,:])):
                raise ValueError("*** Some upper limits are NAN")
            # use 'good' to exclude plotting errorbars for upper limits
            pl.errorbar(eupper.value[good],
                        np.log10(nupperoverg_tofit)[good],
                        yerr=yerr[:,good],
                        linestyle='none',
                        linewidth=0.5,
                        color='k',
                        marker='.', zorder=-5,
                        markersize=2)
        xax = np.array([0, eupper.max().value])
        line = (xax*result.slope.value +
                result.intercept.value)
        pl.plot(xax, np.log10(np.exp(line)), '--',
                color=color,
                alpha=0.6,
                linewidth=1.0,
                label='{2}$T={0:0.1f}$ $\log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value), label))
        pl.ylabel("log N$_u$ (cm$^{-2}$)")
        pl.xlabel("E$_u$ (K)")

        if (uplims is not None) and ((errors is None) or replace_errors_with_uplims):
            # if errors are specified, their errorbars will be better
            # representations of what's actually being fit
            pl.plot(eupper, np.log10(uplims), marker='_', alpha=0.5,
                    linestyle='none', color='k')

    return Ntot, tex, result.slope, result.intercept

def get_deg_(frqs):
    def get_deg(freq):
        closest = np.nanargmin(np.abs(freq-frqs))
        assert frqs[closest]-freq < 1*u.MHz
        return degeneracies[closest]
    return get_deg

def get_Aul_(frqs):
    def get_Aul(freq):
        closest = np.nanargmin(np.abs(freq-frqs))
        assert frqs[closest]-freq < 1*u.MHz
        return einsteinAij[closest]
    return get_Aul

if __name__ == "__main__":

    all_lines = {**lines.disk_lines, **lines.absorbers}

    frequencies = u.Quantity([all_lines[x] for x in all_lines
                              if 'Cl' in x])


    kcl35 = Vamdc.query_molecule('KCl-35')
    rt35 = kcl35.data['RadiativeTransitions']
    frqs = u.Quantity([(float(rt35[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                                u.spectral())
                       for key in rt35])



    upperStateRefs = [rt35[key].UpperStateRef for key in rt35]
    degeneracies = [int(kcl35.data['States'][upperStateRef].TotalStatisticalWeight)
                    for upperStateRef in upperStateRefs]
    einsteinAij = u.Quantity([float(rt35[key].TransitionProbabilityA) for key in rt35], 1/u.s)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    kcl35mask = np.array([(not hasattr(row['Species'], 'mask')) and
                        ('KCl' == row['Species'][:3]) for row in tbl])
    kcl35tbl = tbl[kcl35mask]
    kcl35freqs = u.Quantity(kcl35tbl['Frequency'], u.GHz)
    kkms_kcl35 = (2*np.pi*kcl35tbl['Fitted Width']**2)**0.5 * kcl35tbl['Fitted Amplitude K']
    #kkms_kcl35 = (2*np.pi*(4)**2)**0.5 * kcl35tbl['Fitted Amplitude K']
    Aul = u.Quantity(list(map(get_Aul_(frqs), kcl35freqs)), u.Hz)
    deg = u.Quantity(list(map(get_deg_(frqs), kcl35freqs)))
    kcl35tbl.add_column(table.Column(name='Aul', data=Aul))
    kcl35tbl.add_column(table.Column(name='Degeneracy', data=deg))
    kcl35_nu = nupper_of_kkms(kkms=kkms_kcl35,
                            freq=kcl35freqs,
                            Aul=Aul,
                            degeneracies=deg)
    kcl35tbl.add_column(table.Column(name='N_U', data=kcl35_nu))

    v0 = np.array(['v=0' in row['Species'] for row in kcl35tbl])
    v1 = np.array(['v=1' in row['Species'] for row in kcl35tbl])
    v2 = np.array(['v=2' in row['Species'] for row in kcl35tbl])

    pl.figure(1).clf()
    print("KCl")
    tex0 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v0], u.K), kcl35_nu[v0], plot=True,
                   verbose=True, molecule=kcl35, marker='o', color='r', label='v=0 ')
    tex1 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v1], u.K), kcl35_nu[v1], plot=True,
                   verbose=True, molecule=kcl35, marker='s', color='b', label='v=1 ')
    tex2 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v2], u.K), kcl35_nu[v2], plot=True,
                   verbose=True, molecule=kcl35, marker='^', color='g', label='v=2 ')
    pl.legend(loc='best')
    pl.title("KCl")
    pl.savefig(paths.fpath("KCl_rotational_diagrams.pdf"))




    kcl37 = Vamdc.query_molecule('KCl-37')
    rt37 = kcl37.data['RadiativeTransitions']
    frqs = u.Quantity([(float(rt37[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                                  u.spectral())
                       for key in rt37])



    upperStateRefs = [rt37[key].UpperStateRef for key in rt37]
    degeneracies = [int(kcl37.data['States'][upperStateRef].TotalStatisticalWeight)
                    for upperStateRef in upperStateRefs]
    einsteinAij = u.Quantity([float(rt37[key].TransitionProbabilityA) for key in rt37], 1/u.s)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    v0mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=0' in row['Species'])
                       for row in tbl])
    v1mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=1' in row['Species'])
                       for row in tbl])
    v2mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=2' in row['Species'])
                       for row in tbl])

    kcl37mask = np.array([(not hasattr(row['Species'], 'mask')) and
                          ('K37Cl' == row['Species'][:5])
                          for row in tbl])

    # mask out a bad fit
    bad = (tbl['Species'] == 'K37Clv=1') & (tbl['QNs'] == '29-28')

    kcl37tbl = tbl[kcl37mask & (v0mask | v1mask) & (~bad)]
    kcl37freqs = u.Quantity(kcl37tbl['Frequency'], u.GHz)
    kkms_kcl37 = (2*np.pi*kcl37tbl['Fitted Width']**2)**0.5 * kcl37tbl['Fitted Amplitude K']
    Aul = u.Quantity(list(map(get_Aul_(frqs), kcl37freqs)), u.Hz)
    deg = u.Quantity(list(map(get_deg_(frqs), kcl37freqs)))
    kcl37_nu = nupper_of_kkms(kkms=kkms_kcl37,
                              freq=kcl37freqs,
                              Aul=Aul,
                              degeneracies=deg)

    v0 = np.array(['v=0' in row['Species'] for row in kcl37tbl])
    v1 = np.array(['v=1' in row['Species'] for row in kcl37tbl])

    pl.figure(2).clf()
    print("K 37Cl")
    tex0 = fit_tex(u.Quantity(kcl37tbl['EU_K'][v0], u.K), kcl37_nu[v0], plot=True,
                   verbose=True, molecule=kcl37, marker='o', color='r', label='v=0 ')
    tex1 = fit_tex(u.Quantity(kcl37tbl['EU_K'][v1], u.K), kcl37_nu[v1], plot=True,
                   verbose=True, molecule=kcl37, marker='s', color='b', label='v=1 ')
    pl.legend(loc='best')
    pl.title("K$^{37}$Cl")
    pl.savefig(paths.fpath("KCl37_rotational_diagrams.pdf"))



    nacl = Vamdc.query_molecule('NaCl$')
    rt_nacl = nacl.data['RadiativeTransitions']
    frqs = u.Quantity([(float(rt_nacl[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                                     u.spectral())
                       for key in rt_nacl])



    upperStateRefs = [rt_nacl[key].UpperStateRef for key in rt_nacl]
    degeneracies = [int(nacl.data['States'][upperStateRef].TotalStatisticalWeight)
                    for upperStateRef in upperStateRefs]
    einsteinAij = u.Quantity([float(rt_nacl[key].TransitionProbabilityA) for key in rt_nacl], 1/u.s)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    naclmask = np.array([(not hasattr(row['Species'], 'mask')) and
                         ('NaCl' == row['Species'][:4]) for row in tbl])

    # mask out a bad fit
    # on edge of absorption feature
    bad = (tbl['Species'] == 'NaClv=4') & (tbl['QNs'] == '17-16')


    nacltbl = tbl[naclmask & (~bad)]
    naclfreqs = u.Quantity(nacltbl['Frequency'], u.GHz)
    kkms_nacl = (2*np.pi*nacltbl['Fitted Width']**2)**0.5 * nacltbl['Fitted Amplitude K']
    Aul = u.Quantity(list(map(get_Aul_(frqs), naclfreqs)), u.Hz)
    deg = u.Quantity(list(map(get_deg_(frqs), naclfreqs)))
    nacl_nu = nupper_of_kkms(kkms=kkms_nacl,
                             freq=naclfreqs,
                             Aul=Aul,
                             degeneracies=deg)


    #v0 = np.array(['v=0' in row['Species'] for row in nacltbl])
    v1 = np.array(['v=1' in row['Species'] for row in nacltbl])
    v2 = np.array(['v=2' in row['Species'] for row in nacltbl])
    v3 = np.array(['v=3' in row['Species'] for row in nacltbl])
    v4 = np.array(['v=4' in row['Species'] for row in nacltbl])

    pl.figure(3).clf()
    print("NaCl")
    #tex0 = fit_tex(u.Quantity(nacltbl['EU_K'][v0], u.K), nacl_nu[v0], plot=True,
    #               verbose=True, molecule=nacl, marker='o', color='r')
    tex1 = fit_tex(u.Quantity(nacltbl['EU_K'][v1], u.K), nacl_nu[v1], plot=True,
                   verbose=True, molecule=nacl, marker='s', color='b', label='v=1 ')
    tex2 = fit_tex(u.Quantity(nacltbl['EU_K'][v2], u.K), nacl_nu[v2], plot=True,
                   verbose=True, molecule=nacl, marker='^', color='g', label='v=2 ')
    tex3 = fit_tex(u.Quantity(nacltbl['EU_K'][v3], u.K), nacl_nu[v3], plot=True,
                   verbose=True, molecule=nacl, marker='o', color='r', label='v=3 ')
    tex4 = fit_tex(u.Quantity(nacltbl['EU_K'][v4], u.K), nacl_nu[v4], plot=True,
                   verbose=True, molecule=nacl, marker='d', color='orange', label='v=4 ')
    pl.legend(loc='best')
    pl.axis([300,1600,9.0,13])
    pl.title("NaCl")
    pl.savefig(paths.fpath("NaCl_rotational_diagrams.pdf"))








    nacl37 = Vamdc.query_molecule('NaCl-37')
    rt_nacl37 = nacl37.data['RadiativeTransitions']
    frqs = u.Quantity([(float(rt_nacl37[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                                     u.spectral())
                       for key in rt_nacl37])



    upperStateRefs = [rt_nacl37[key].UpperStateRef for key in rt_nacl37]
    degeneracies = [int(nacl37.data['States'][upperStateRef].TotalStatisticalWeight)
                    for upperStateRef in upperStateRefs]
    einsteinAij = u.Quantity([float(rt_nacl37[key].TransitionProbabilityA) for key in rt_nacl37], 1/u.s)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    nacl37mask = np.array([(not hasattr(row['Species'], 'mask')) and
                         ('Na37Cl' == row['Species'][:6]) for row in tbl])
    nacl37tbl = tbl[nacl37mask & (v0mask | v1mask)]
    nacl37freqs = u.Quantity(nacl37tbl['Frequency'], u.GHz)
    kkms_nacl37 = (2*np.pi*nacl37tbl['Fitted Width']**2)**0.5 * nacl37tbl['Fitted Amplitude K']
    Aul = u.Quantity(list(map(get_Aul_(frqs), nacl37freqs)), u.Hz)
    deg = u.Quantity(list(map(get_deg_(frqs), nacl37freqs)))
    nacl37_nu = nupper_of_kkms(kkms=kkms_nacl37,
                               freq=nacl37freqs,
                               Aul=Aul,
                               degeneracies=deg)


    v0 = np.array(['v=0' in row['Species'] for row in nacl37tbl])
    v1 = np.array(['v=1' in row['Species'] for row in nacl37tbl])
    #v2 = np.array(['v=2' in row['Species'] for row in nacl37tbl])
    #v3 = np.array(['v=3' in row['Species'] for row in nacl37tbl])
    #v4 = np.array(['v=4' in row['Species'] for row in nacl37tbl])

    pl.figure(4).clf()
    print("Na37Cl")
    tex0 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v0], u.K), nacl37_nu[v0], plot=True,
                   verbose=True, molecule=nacl37, marker='o', color='r')
    tex1 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v1], u.K), nacl37_nu[v1], plot=True,
                   verbose=True, molecule=nacl37, marker='s', color='b', label='v=1 ')
    #tex2 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v2], u.K), nacl37_nu[v2], plot=True,
    #               verbose=True, molecule=nacl37, marker='^', color='g', label='v=2 ')
    #tex3 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v3], u.K), nacl37_nu[v3], plot=True,
    #               verbose=True, molecule=nacl37, marker='o', color='r', label='v=3 ')
    #tex4 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v4], u.K), nacl37_nu[v4], plot=True,
    #               verbose=True, molecule=nacl37, marker='d', color='orange', label='v=4 ')
    pl.legend(loc='best')
    pl.title("Na$^{37}$Cl")
    pl.axis([-50,650,10.5,12.15])
    pl.savefig(paths.fpath("NaCl37_rotational_diagrams.pdf"))






    k41cl = Vamdc.query_molecule('K-41-Cl')
    rt41 = k41cl.data['RadiativeTransitions']
    frqs = u.Quantity([(float(rt41[key].FrequencyValue)*u.MHz).to(u.GHz,
                                                                  u.spectral())
                       for key in rt41])



    upperStateRefs = [rt41[key].UpperStateRef for key in rt41]
    degeneracies = [int(k41cl.data['States'][upperStateRef].TotalStatisticalWeight)
                    for upperStateRef in upperStateRefs]
    einsteinAij = u.Quantity([float(rt41[key].TransitionProbabilityA) for key in rt41], 1/u.s)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    k41clmask = np.array([(not hasattr(row['Species'], 'mask')) and
                          ('41KCl' == row['Species'][:5])
                          for row in tbl])

    bad = (((tbl['Species'] == '41KClv=2') & (tbl['QNs'] == '29-28')) | # absorption
           ((tbl['Species'] == '41KClv=0') & (tbl['QNs'] == '31-30'))) # wing of NaCl

    k41cltbl = tbl[k41clmask & (v0mask | v1mask) & (~bad)]
    k41clfreqs = u.Quantity(k41cltbl['Frequency'], u.GHz)
    kkms_k41cl = (2*np.pi*k41cltbl['Fitted Width']**2)**0.5 * k41cltbl['Fitted Amplitude K']
    Aul = u.Quantity(list(map(get_Aul_(frqs), k41clfreqs)), u.Hz)
    deg = u.Quantity(list(map(get_deg_(frqs), k41clfreqs)))
    k41cl_nu = nupper_of_kkms(kkms=kkms_k41cl,
                              freq=k41clfreqs,
                              Aul=Aul,
                              degeneracies=deg)

    v0 = np.array(['v=0' in row['Species'] for row in k41cltbl])
    v1 = np.array(['v=1' in row['Species'] for row in k41cltbl])

    pl.figure(5).clf()
    print("41KCl")
    tex0 = fit_tex(u.Quantity(k41cltbl['EU_K'][v0], u.K), k41cl_nu[v0], plot=True,
                   verbose=True, molecule=k41cl, marker='o', color='r', label='v=0 ')
    tex1 = fit_tex(u.Quantity(k41cltbl['EU_K'][v1], u.K), k41cl_nu[v1], plot=True,
                   verbose=True, molecule=k41cl, marker='s', color='b', label='v=1 ')
    pl.legend(loc='best')
    pl.title("$^{41}$KCl")
    pl.savefig(paths.fpath("K41Cl_rotational_diagrams.pdf"))

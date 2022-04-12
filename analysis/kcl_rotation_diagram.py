import paths
import numpy as np
from astropy import units as u
from astropy import constants
from astropy import log
from astropy import table
from astropy import modeling
from astropy.modeling.models import custom_model

import dust_emissivity


import pylab as pl

import lines
import salt_tables
from lte_modeling_tools import (rovib_lte_model_generator,
                                simple_lte_model_generator, kkms_of_nupper,
                                nupper_of_kkms)
from pyspeckit.spectrum.models import lte_molecule as plm
from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters

vib_constants = {'KCl': (281*u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()),
                 'K37Cl': (275.87497*u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()),
                 '41KCl': (276.639613*u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()),
                 'NaCl': (366*u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()),
                 # couldn't find numbers; assume it's reduced by mass (which
                 # isn't a good assumption)
                 'Na37Cl': (35+23)/(37+23)*(366*u.cm**-1).to(u.eV, u.spectral()).to(u.K, u.temperature_energy()),
                }

def get_vib_energies(moltbl):
    return {v: moltbl['E_L'][(moltbl['vu'] == v) & (moltbl['Jl'] == 0)].quantity[0].value
            for v in np.unique(moltbl['vu'])}

def get_rot_energies(moltbl):
    def matchJu(ju):
        rslt = moltbl['E_U'][(moltbl['vu'] == 0) & (moltbl['Ju'] == ju)].quantity
        if len(rslt) == 1:
            return rslt
    return {ju: matchJu(ju)[0].value
            for ju in np.unique(moltbl['Ju'])
            if matchJu(ju)
           }


def fit_multi_tex(eupper, nupperoverg, vstate, jstate, vibenergies,
                  rotenergies, verbose=False, plot=False, uplims=None,
                  errors=None,
                  min_nupper=1,
                  replace_errors_with_uplims=False,
                  #molecule=None,
                  partition_func=None,
                  colors='rgbcmyk',
                  molname=None,
                  collims=(np.log(1e9), np.log(1e17)),
                  rottemlims=(10,300),
                  vibtemlims=(500,16000),
                  marker='o', max_uplims='half'):
    """
    Fit the Boltzmann diagram with a vibrational and a rotational temperature

    Parameters
    ----------
    max_uplims: str or number
        The maximum number of upper limits before the fit is ignored completely
        and instead zeros are returned
    """

    nupperoverg_tofit = nupperoverg.copy().to(u.cm**-2).value
    if errors is not None:
        errors = errors.to(u.cm**-2).value

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

    model = rovib_lte_model_generator(vibenergies, rotenergies)

    #print('vibenergies:',vibenergies, '\nrotenergies:', rotenergies)

    fitter = modeling.fitting.LevMarLSQFitter()
    #fitter = modeling.fitting.LinearLSQFitter()

    model_to_fit = model
    model_to_fit.logcolumn.bounds = collims
    model_to_fit.rottem.bounds = rottemlims
    model_to_fit.vibtem.bounds = vibtemlims
    #model_to_fit.rottem.fixed = True

    # check that the model works
    #first_check = model_to_fit(eupper[good].to(u.K).value)
    #print("Result of first check: ",first_check)

    result = fitter(model_to_fit,
                    jstate[good].astype('int'),
                    vstate[good].astype('int'),
                    np.log(nupperoverg_tofit[good]),
                    weights=np.log(weights[good]))
    print(result)
    tex = result.rottem #u.Quantity(-1./result.slope, u.K)
    tvib = result.vibtem

    #partition_func = specmodel.calculate_partitionfunction(molecule.data['States'],
    #                                                       temperature=tex.value)
    #assert len(partition_func) == 1
    #Q_rot = tuple(partition_func.values())[0]
    #print("Q_rot:", Q_rot)
    Q_rot = partition_func(tex.value)

    Ntot = np.exp(result.logcolumn + np.log(Q_rot)) * u.cm**-2

    if verbose:
        print(("Tex={0}, Ntot={1}, log(Ntot)={4}, Q_rot={2}, "
               "nuplim={3}".format(tex, Ntot, Q_rot, upperlim_mask.sum(),
                                   np.log10(Ntot.value),
                                  )))

    if plot:
        import pylab as pl
        for vib,color in zip(np.arange(vstate.max()+1), colors):
            mask = vstate == vib
            if mask.sum() == 0:
                continue
            L, = pl.plot(eupper[mask], np.log10(nupperoverg_tofit[mask]),
                         marker=marker,
                         color=color, markeredgecolor='none', alpha=0.5,
                         linestyle='none',
                         #markersize=2,
                        )
            if uplims is not None:
                L, = pl.plot(eupper[upperlim_mask & mask],
                             np.log10(uplims)[upperlim_mask & mask], 'bv', alpha=0.2,
                             markersize=2)
                #L, = pl.plot(eupper[upperlim_mask],
                #             np.log10(nupperoverg)[upperlim_mask], 'bv', alpha=0.2)
            if errors is not None:
                yerr = np.array([np.log10(nupperoverg_tofit)-np.log10(nupperoverg_tofit-errors),
                                 np.log10(nupperoverg_tofit+errors)-np.log10(nupperoverg_tofit)])
                # if lower limit is nan, set to zero
                yerr[0,:] = np.nan_to_num(yerr[0,:])
                if np.any(np.isnan(yerr[1,:])):
                    print(ValueError("*** Some upper limits are NAN"))
                # use 'good' to exclude plotting errorbars for upper limits
                pl.errorbar(eupper.value[good & mask],
                            np.log10(nupperoverg_tofit)[good & mask],
                            yerr=yerr[:,good & mask],
                            linestyle='none',
                            linewidth=0.5,
                            color=L.get_color(),
                            marker='.', zorder=-5,
                            markersize=2)
            #xax = np.array([0, eupper.max().value+500])

            jstates = np.arange(1,np.max(list(rotenergies.keys()))+1)
            vstates = np.ones(np.max(list(rotenergies.keys())))*vib
            line = result(jstates, vstates)
            xax = np.array([rotenergies[ju] + vibenergies[vib]
                            for ju in jstates])

            for marker, mask in (('.', (jstates>-1) & (jstates <= 10)),
                                 ('s', (jstates>10) & (jstates <= 20)),
                                 ('o', (jstates>20) & (jstates <= 30)),
                                 ('D', (jstates>30) & (jstates <= 40)),
                                 ('p', (jstates>40) & (jstates <= 50))):
                pl.plot(xax[mask],
                        np.log10(np.exp(line))[mask],
                        marker=marker,
                        linestyle='none',
                        color=color,
                        alpha=0.3,
                        linewidth=1.0,
                        markersize=4,
                       )
                    #label='v={0}'.format(vib),
                    #label=('v={2}$T_R={0:0.1f}$ $T_v={3:0.1f}$ $\log(N)={1:0.1f}$'
                    #       .format(tex.value, np.log10(Ntot.value), vib,
                    #               tvib.value
                    #              ))

        for jj in np.unique(jstate):
            #jstates = np.arange(1,np.max(list(rotenergies.keys()))+1)
            #vstates = np.ones(np.max(list(rotenergies.keys())))*vib
            vstates = np.arange(9)
            jstates = np.ones(9)*jj
            line = result(jstates, vstates)
            xax = np.array([rotenergies[jj] + vibenergies[vv]
                            for vv in vstates])

            pl.plot(xax, np.log10(np.exp(line)), '--',
                    color='k',
                    alpha=0.1,
                    linewidth=1.0,
                   )
                    #label="Ju={0}".format(jj),
                    #label=('v={2}$T_R={0:0.1f}$ $T_v={3:0.1f}$ $\log(N)={1:0.1f}$'
                    #       .format(tex.value, np.log10(Ntot.value), vib,
                    #               tvib.value
                    #              ))

        for eu, nu, jj, vv in zip(eupper, nupperoverg_tofit, jstate, vstate):
            model_nu = result(jj, vv)
            #print("POINT MATCHING: ",eu, model_nu, np.log(nu), jj, vv)
            pl.plot(u.Quantity([eu, eu]),
                    np.log10(np.exp([np.log(nu), model_nu])),
                    linestyle='-', alpha=0.3,
                    color='k',
                    linewidth=0.5)


        pl.ylabel("log N$_u$/g (cm$^{-2}$)")
        pl.xlabel("E$_u$ (K)")

        pl.plot([], [], linestyle='-', color='k',
                label=('$T_R={0:0.1f}$ $T_v={2:0.1f}$ $\log(N)={1:0.1f}$'
                       .format(tex.value, np.log10(Ntot.value), tvib.value)))

        if (uplims is not None) and ((errors is None) or replace_errors_with_uplims):
            # if errors are specified, their errorbars will be better
            # representations of what's actually being fit
            pl.plot(eupper, np.log10(uplims), marker='_', alpha=0.5,
                    linestyle='none', color='k')

    return Ntot, tex, result.rottem, result.vibtem, result.logcolumn

def fit_tex(eupper, nupperoverg, verbose=False, plot=False, uplims=None,
            errors=None, min_nupper=1,
            replace_errors_with_uplims=False,
            #molecule=None,
            partition_func=None,
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

    nupperoverg_tofit = nupperoverg.copy().to(u.cm**-2).value
    if errors is not None:
        errors = errors.to(u.cm**-2).value

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

    #model = modeling.models.Linear1D()
    model = simple_lte_model_generator()
    fitter = modeling.fitting.LevMarLSQFitter()
    #fitter = modeling.fitting.LinearLSQFitter()
    #model.slope.bounds = [0, -1000]
    if good.sum() >= 2:
        result = fitter(model,
                        eupper[good].to(u.K).value,
                        np.log(nupperoverg_tofit[good]),
                        weights=np.log(weights[good]))
        #tex = u.Quantity(-1./result.slope, u.K)
        tex = u.Quantity(result.tem, u.K)
        if tex < 0*u.K:
            tex = 100 * u.K

        #partition_func = specmodel.calculate_partitionfunction(molecule.data['States'],
        #                                                       temperature=tex.value)
        #assert len(partition_func) == 1
        #Q_rot = tuple(partition_func.values())[0]
        Q_rot = partition_func(tex.value)

        #Ntot = np.exp(result.intercept + np.log(Q_rot)) * u.cm**-2
        Ntot = np.exp(result.logcolumn) * Q_rot * u.cm**-2
    else:
        Ntot = 1e10*u.cm**-2
        tex = 100*u.K
        Q_rot = 1
        result = model

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
                #raise ValueError("*** Some upper limits are NAN")
                print(ValueError("*** Some upper limits are NAN"))
            # use 'good' to exclude plotting errorbars for upper limits
            pl.errorbar(eupper.value[good],
                        np.log10(nupperoverg_tofit)[good],
                        yerr=yerr[:,good],
                        linestyle='none',
                        linewidth=0.5,
                        color='k',
                        marker='.', zorder=-5,
                        markersize=2)
        xax = np.array([0, eupper.max().value+500])
        line = (xax*(-1/result.tem.value) +
                result.logcolumn.value)
        pl.plot(xax, np.log10(np.exp(line)), '--',
                color=color,
                alpha=0.6,
                linewidth=1.0,
                label='{2}$T={0:0.1f}$ $\log(N)={1:0.1f}$'.format(tex, np.log10(Ntot.value), label))
        pl.ylabel("log N$_u$/g (cm$^{-2}$)")
        pl.xlabel("E$_u$ (K)")

        if (uplims is not None) and ((errors is None) or replace_errors_with_uplims):
            # if errors are specified, their errorbars will be better
            # representations of what's actually being fit
            pl.plot(eupper, np.log10(uplims), marker='_', alpha=0.5,
                    linestyle='none', color='k')

    return Ntot, tex, result.tem, result.logcolumn

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


    # from astroquery.vamdc import Vamdc
    # from vamdclib import specmodel

    # kcl35 = Vamdc.query_molecule('KCl-35')
    # rt35 = kcl35.data['RadiativeTransitions']
    # frqs = u.Quantity([(float(rt35[key].FrequencyValue)*u.MHz).to(u.GHz,
    #                                                               u.spectral())
    #                    for key in rt35])


    # upperStateRefs = [rt35[key].UpperStateRef for key in rt35]
    # degeneracies = [int(kcl35.data['States'][upperStateRef].TotalStatisticalWeight)
    #                 for upperStateRef in upperStateRefs]
    # einsteinAij = u.Quantity([float(rt35[key].TransitionProbabilityA) for key in rt35], 1/u.s)


    frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('KCl', fmin=85*u.GHz, fmax=360*u.GHz)


    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    kcl35mask = np.array([(not hasattr(row['Species'], 'mask')) and
                         ('KCl' == row['Species'][:3] or
                          '39K-35Cl' in row['Species']) for row in tbl])

    bad = ((tbl['Line Name'] == 'KClv=4_13-12') |
           (tbl['Line Name'] == 'KClv=6_31-30') |
           (tbl['Line Name'] == 'KClv=8_27-26'))

    print("KCl: {0} in-band, {1} detected".format(kcl35mask.sum(),
                                                  (kcl35mask & (~bad)).sum()))

    kcl35tbl = tbl[kcl35mask]
    kcl35freqs = u.Quantity(kcl35tbl['Frequency'], u.GHz)
    kkms_kcl35 = (2*np.pi*kcl35tbl['Fitted Width']**2)**0.5 * kcl35tbl['Fitted Amplitude K']
    ekkms_kcl35 = (2*np.pi)**0.5 * (kcl35tbl['Fitted Width error']**2 *
                                    kcl35tbl['Fitted Amplitude K']**2 +
                                    kcl35tbl['Fitted Width']**2 *
                                    kcl35tbl['Fitted Amplitude error K']**2)**0.5
    #kkms_kcl35 = (2*np.pi*(4)**2)**0.5 * kcl35tbl['Fitted Amplitude K']
    #Aul = u.Quantity(list(map(get_Aul_(frqs), kcl35freqs)), u.Hz)
    #deg = u.Quantity(list(map(get_deg_(frqs), kcl35freqs)))
    #kcl35tbl.add_column(table.Column(name='Aul', data=Aul))
    #kcl35tbl.add_column(table.Column(name='Degeneracy', data=deg))
    Aul = kcl35tbl['Aij']
    deg = kcl35tbl['deg']
    kcl35_nug = nupper_of_kkms(kkms=kkms_kcl35,
                               freq=kcl35freqs,
                               Aul=Aul,
                               #degeneracies=deg,
                               )
    kcl35_nu = kcl35_nug / deg
    ekcl35_nug = nupper_of_kkms(kkms=ekkms_kcl35,
                                freq=kcl35freqs,
                                Aul=Aul,
                                #degeneracies=deg
                               )
    ekcl35_nu = ekcl35_nug / deg
    kcl35tbl.add_column(table.Column(name='N_U', data=kcl35_nu,
                                     ))

    v0 = np.array(['v=0' in row['Species'] for row in kcl35tbl])
    v1 = np.array(['v=1' in row['Species'] for row in kcl35tbl])
    v2 = np.array(['v=2' in row['Species'] for row in kcl35tbl])
    v3 = np.array(['v=3' in row['Species'] for row in kcl35tbl])
    v4 = np.array(['v=4' in row['Species'] for row in kcl35tbl])
    v5 = np.array(['v=5' in row['Species'] for row in kcl35tbl])
    v6 = np.array(['v=6' in row['Species'] for row in kcl35tbl])
    v7 = np.array(['v=7' in row['Species'] for row in kcl35tbl])
    v8 = np.array(['v=8' in row['Species'] for row in kcl35tbl])
    j13 = np.array(['J=13' in row['Species'] for row in kcl35tbl])
    j29 = np.array(['J=29' in row['Species'] for row in kcl35tbl])
    j31 = np.array(['J=31' in row['Species'] for row in kcl35tbl])

    pl.figure(1).clf()
    print("KCl")

    optdepth = (dust_emissivity.dust.kappa(kcl35freqs, beta=2) * 120*u.g/u.cm**2).decompose().value
    kcl35_nug_dustcorr = nupper_of_kkms(kkms=kkms_kcl35 * np.exp(optdepth),
                                        freq=kcl35freqs,
                                        Aul=Aul,
                                        #degeneracies=deg,
                                        )
    kcl35_nu_dustcorr = kcl35_nug_dustcorr / deg

    tex0 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v0], u.K), kcl35_nu[v0],
                   errors=ekcl35_nu[v0], plot=True, verbose=True,
                   #partition_func=partfunc,
                   partition_func=partfunc,
                   marker='o', color='r', label='v=0 ')
    #tex0_dustc = fit_tex(u.Quantity(kcl35tbl['EU_K'][v0], u.K), kcl35_nu_dustcorr[v0],
    #                     errors=ekcl35_nu[v0], plot=True, verbose=True,
    #                     partition_func=partfunc, marker='o', color='b', label='v=0 ')
    #tex0_dustc_nob3 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v0][1:], u.K), kcl35_nu_dustcorr[v0][1:],
    #                          errors=ekcl35_nu[v0][1:], plot=True, verbose=True,
    #                     partition_func=partfunc, marker='o', color='g', label='v=0 no b3')
    #tex1 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v1], u.K), kcl35_nu[v1],
    #               errors=ekcl35_nu[v1], plot=True, verbose=True,
    #               partition_func=partfunc, marker='s', color='b', label='v=1 ')
    #tex2 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v2], u.K), kcl35_nu[v2],
    #               errors=ekcl35_nu[v2], plot=True, verbose=True,
    #               partition_func=partfunc, marker='^', color='g', label='v=2 ')
    #tex3 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v3], u.K), kcl35_nu[v3],
    #               errors=ekcl35_nu[v3], plot=True, verbose=True,
    #               partition_func=partfunc, marker='v', color='orange', label='v=3 ')
    #tex4 = fit_tex(u.Quantity(kcl35tbl['EU_K'][v4], u.K), kcl35_nu[v4],
    #               errors=ekcl35_nu[v4], plot=True, verbose=True,
    #               partition_func=partfunc, marker='d', color='m', label='v=4 ')
    texJ13 = fit_tex(u.Quantity(kcl35tbl['EU_K'][j13], u.K), kcl35_nu[j13],
                     errors=ekcl35_nu[j13], plot=True, verbose=True,
                     partition_func=partfunc, marker='h', color='k', label='J=13 ')
    #texJ13_dustc = fit_tex(u.Quantity(kcl35tbl['EU_K'][j13], u.K), kcl35_nu_dustcorr[j13],
    #                 errors=ekcl35_nu[j13], plot=True, verbose=True,
    #                 partition_func=partfunc, marker='h', color='k', label='J=13 dust')
    texJ29 = fit_tex(u.Quantity(kcl35tbl['EU_K'][j29], u.K), kcl35_nu[j29],
                     errors=ekcl35_nu[j29], plot=True, verbose=True,
                     partition_func=partfunc, marker='o', color='b', label='J=29 ')
    texJ31 = fit_tex(u.Quantity(kcl35tbl['EU_K'][j31], u.K), kcl35_nu[j31],
                     errors=ekcl35_nu[j31], plot=True, verbose=True,
                     partition_func=partfunc, marker='P', color='navy', label='J=31 ')

    pl.legend(loc='lower right')
    pl.title("KCl")
    pl.axis((0, 2500, 7.75, 10.5))
    pl.savefig(paths.fpath("KCl_rotational_diagrams.pdf"))


    pl.figure(1).clf()
    vstate = 0*v0 + 1*v1 + 2*v2 + 3*v3 + 4*v4 + 5*v5 + 6*v6 + 7*v7 + 8*v8

    print(fit_multi_tex(eupper=u.Quantity(kcl35tbl['EU_K'], u.K),
                        nupperoverg=kcl35_nu,
                        vstate=vstate,
                        jstate=np.array(kcl35tbl['J$_u$'], dtype='int'),
                        vibenergies=get_vib_energies(salt_tables.KCl),
                        rotenergies=get_rot_energies(salt_tables.KCl),
                        errors=ekcl35_nu,
                        rottemlims=(30,150),
                        vibtemlims=(2500,5000),
                        collims=(np.log(1e8), np.log(1e16)),
                        plot=True, verbose=True, partition_func=partfunc, marker='o',
                        molname='KCl',
                        colors=('r','g','b','orange','m','c','darkred','darkgreen',
                                'purple'), )
         )
    pl.legend(loc='lower right')
    pl.title("KCl")
    pl.axis((0, 4000, 8.0, 10.5))

    pl.savefig(paths.fpath("KCl_rotational-vibrational_fit_diagrams.pdf"))




    #kcl37 = Vamdc.query_molecule('KCl-37')
    #rt37 = kcl37.data['RadiativeTransitions']
    #frqs = u.Quantity([(float(rt37[key].FrequencyValue)*u.MHz).to(u.GHz,
    #                                                              u.spectral())
    #                   for key in rt37])



    #upperStateRefs = [rt37[key].UpperStateRef for key in rt37]
    #degeneracies = [int(kcl37.data['States'][upperStateRef].TotalStatisticalWeight)
    #                for upperStateRef in upperStateRefs]
    #einsteinAij = u.Quantity([float(rt37[key].TransitionProbabilityA) for key in rt37], 1/u.s)

    frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('KCl-37, v=0-15', catalog='CDMS', fmin=85*u.GHz, fmax=360*u.GHz)

    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    v0mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=0' in row['Species'])
                       for row in tbl])
    v1mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=1' in row['Species'])
                       for row in tbl])
    v2mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=2' in row['Species'])
                       for row in tbl])
    v3mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=3' in row['Species'])
                       for row in tbl])
    v4mask = np.array([(not hasattr(row['Species'], 'mask')) and ('v=4' in row['Species'])
                       for row in tbl])

    kcl37mask = np.array([(not hasattr(row['Species'], 'mask')) and
                          ('K37Cl' == row['Species'][:5] or
                           '39K-37Cl' in row['Species'])
                          for row in tbl])

    # mask out a bad fit
    bad = ((tbl['Line Name'] == 'K37Clv=1_29-28') | 
           (tbl['Line Name'] == '41KClv=0 46-45'))

    print("K37Cl: {0} in-band, {1} detected".format(kcl37mask.sum(),
                                                    (kcl37mask & (~bad)).sum()))

    kcl37tbl = tbl[kcl37mask & (~bad)]
    kcl37freqs = u.Quantity(kcl37tbl['Frequency'], u.GHz)
    kkms_kcl37 = (2*np.pi*kcl37tbl['Fitted Width']**2)**0.5 * kcl37tbl['Fitted Amplitude K']
    ekkms_kcl37 = (2*np.pi)**0.5 * (kcl37tbl['Fitted Width error']**2 *
                                    kcl37tbl['Fitted Amplitude K']**2 +
                                    kcl37tbl['Fitted Width']**2 *
                                    kcl37tbl['Fitted Amplitude error K']**2)**0.5
    #Aul = u.Quantity(list(map(get_Aul_(frqs), kcl37freqs)), u.Hz)
    #deg = u.Quantity(list(map(get_deg_(frqs), kcl37freqs)))
    Aul = kcl37tbl['Aij']
    deg = kcl37tbl['deg']
    kcl37_nu = nupper_of_kkms(kkms=kkms_kcl37,
                              freq=kcl37freqs,
                              Aul=Aul)/deg
    ekcl37_nu = nupper_of_kkms(kkms=ekkms_kcl37,
                               freq=kcl37freqs,
                               Aul=Aul)/deg

    v0 = np.array(['v=0' in row['Species'] for row in kcl37tbl])
    v1 = np.array(['v=1' in row['Species'] for row in kcl37tbl])

    pl.figure(2).clf()
    print("K 37Cl")
    tex0 = fit_tex(u.Quantity(kcl37tbl['EU_K'][v0], u.K), kcl37_nu[v0],
                   errors=ekcl37_nu[v0], plot=True, verbose=True,
                   partition_func=partfunc, marker='o', color='r', label='v=0 ')
    tex1 = fit_tex(u.Quantity(kcl37tbl['EU_K'][v1], u.K), kcl37_nu[v1],
                   errors=ekcl37_nu[v1], plot=True, verbose=True,
                   partition_func=partfunc, marker='s', color='b', label='v=1 ')
    pl.legend(loc='lower right')
    pl.title("K$^{37}$Cl")
    pl.savefig(paths.fpath("KCl37_rotational_diagrams.pdf"))

    pl.figure(2).clf()
    vstate = 0*v0 + 1*v1

    print(fit_multi_tex(eupper=u.Quantity(kcl37tbl['EU_K'], u.K),
                        nupperoverg=kcl37_nu,
                        vstate=vstate,
                        jstate=np.array(kcl37tbl['J$_u$'], dtype='int'),
                        vibenergies=get_vib_energies(salt_tables.K37Cl),
                        rotenergies=get_rot_energies(salt_tables.K37Cl),
                        errors=ekcl37_nu,
                        plot=True, verbose=True, partition_func=partfunc, marker='o',
                        molname='K37Cl',
                        colors=('r','g','b','orange','m','c'), )
         )
    pl.legend(loc='lower right')
    pl.title("K$^{37}$Cl")

    pl.savefig(paths.fpath("KCl37_rotational-vibrational_fit_diagrams.pdf"))


    #nacl = Vamdc.query_molecule('NaCl$')
    #rt_nacl = nacl.data['RadiativeTransitions']
    #frqs = u.Quantity([(float(rt_nacl[key].FrequencyValue)*u.MHz).to(u.GHz,
    #                                                                 u.spectral())
    #                   for key in rt_nacl])



    #upperStateRefs = [rt_nacl[key].UpperStateRef for key in rt_nacl]
    #degeneracies = [int(nacl.data['States'][upperStateRef].TotalStatisticalWeight)
    #                for upperStateRef in upperStateRefs]
    #einsteinAij = u.Quantity([float(rt_nacl[key].TransitionProbabilityA) for key in rt_nacl], 1/u.s)

    frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('NaCl, v=0-15', catalog='CDMS', fmin=85*u.GHz, fmax=360*u.GHz)

    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    naclmask = np.array([(not hasattr(row['Species'], 'mask')) and
                         ('NaCl' == row['Species'][:4] or
                         '23Na-35Cl' in row['Species'])
                         for row in tbl])

    # mask out a bad fit
    # on edge of absorption feature
    bad = (tbl['Line Name'] == 'NaClv=4_17-16')

    print("NaCl: {0} in-band, {1} detected".format(naclmask.sum(),
                                                   (naclmask & (~bad)).sum()))

    nacltbl = tbl[naclmask & (~bad)]
    naclfreqs = u.Quantity(nacltbl['Frequency'], u.GHz)
    kkms_nacl = (2*np.pi*nacltbl['Fitted Width']**2)**0.5 * nacltbl['Fitted Amplitude K']
    ekkms_nacl = (2*np.pi)**0.5 * (nacltbl['Fitted Width error']**2 *
                                   nacltbl['Fitted Amplitude K']**2 +
                                   nacltbl['Fitted Width']**2 *
                                   nacltbl['Fitted Amplitude error K']**2)**0.5
    #Aul = u.Quantity(list(map(get_Aul_(frqs), naclfreqs)), u.Hz)
    #deg = u.Quantity(list(map(get_deg_(frqs), naclfreqs)))
    Aul = nacltbl['Aij']
    deg = nacltbl['deg']
    nacl_nu = nupper_of_kkms(kkms=kkms_nacl,
                             freq=naclfreqs,
                             Aul=Aul,)/deg
    enacl_nu = nupper_of_kkms(kkms=ekkms_nacl,
                              freq=naclfreqs,
                              Aul=Aul)/deg


    v0 = np.array(['v=0' in row['Species'] for row in nacltbl])
    v1 = np.array(['v=1' in row['Species'] for row in nacltbl])
    v2 = np.array(['v=2' in row['Species'] for row in nacltbl])
    v3 = np.array(['v=3' in row['Species'] for row in nacltbl])
    v4 = np.array(['v=4' in row['Species'] for row in nacltbl])
    v5 = np.array(['v=5' in row['Species'] for row in nacltbl])
    v6 = np.array(['v=6' in row['Species'] for row in nacltbl])
    v7 = np.array(['v=7' in row['Species'] for row in nacltbl])
    v8 = np.array(['v=8' in row['Species'] for row in nacltbl])
    j7 = np.array(['J=7' in row['Species'] for row in nacltbl])
    j8 = np.array(['J=8' in row['Species'] for row in nacltbl])

    pl.figure(3).clf()
    print("NaCl")
    #tex0 = fit_tex(u.Quantity(nacltbl['EU_K'][v0], u.K), nacl_nu[v0], plot=True,
    #               verbose=True, partition_func=partfunc, marker='o', color='r')
    tex1 = fit_tex(u.Quantity(nacltbl['EU_K'][v1], u.K), nacl_nu[v1],
                   errors=enacl_nu[v1], plot=True,
                   verbose=True, partition_func=partfunc, marker='s', color='b', label='v=1 ')
    tex2 = fit_tex(u.Quantity(nacltbl['EU_K'][v2], u.K), nacl_nu[v2],
                   errors=enacl_nu[v2], plot=True,
                   verbose=True, partition_func=partfunc, marker='^', color='g', label='v=2 ')
    #tex3 = fit_tex(u.Quantity(nacltbl['EU_K'][v3], u.K), nacl_nu[v3],
    #               errors=enacl_nu[v3], plot=True,
    #               verbose=True, partition_func=partfunc, marker='o', color='r', label='v=3 ')
    #tex4 = fit_tex(u.Quantity(nacltbl['EU_K'][v4], u.K), nacl_nu[v4],
    #               errors=enacl_nu[v4], plot=True,
    #               verbose=True, partition_func=partfunc, marker='d', color='orange', label='v=4 ')
    #tex5 = fit_tex(u.Quantity(nacltbl['EU_K'][v5], u.K), nacl_nu[v5],
    #               errors=enacl_nu[v5], plot=True,
    #               verbose=True, partition_func=partfunc, marker='v', color='m', label='v=5 ')
    #tex6 = fit_tex(u.Quantity(nacltbl['EU_K'][v6], u.K), nacl_nu[v6],
    #               errors=enacl_nu[v6], plot=True,
    #               verbose=True, partition_func=partfunc, marker='v', color='m', label='v=6 ')
    texJ8 = fit_tex(u.Quantity(nacltbl['EU_K'][j8], u.K), nacl_nu[j8],
                    errors=enacl_nu[j8], plot=True, verbose=True,
                    partition_func=partfunc, marker='h', color='k', label='J=8 ')
    texJ7 = fit_tex(u.Quantity(nacltbl['EU_K'][j7], u.K), nacl_nu[j7],
                    errors=enacl_nu[j7], plot=True, verbose=True,
                    partition_func=partfunc, marker='>', color='c', label='J=7 ')
    pl.legend(loc='upper right')
    pl.axis([300,3300,9.0,13])
    pl.title("NaCl")
    pl.savefig(paths.fpath("NaCl_rotational_diagrams.pdf"))


    pl.figure(3).clf()
    vstate = 0*v0 + 1*v1 + 2*v2 + 3*v3 + 4*v4 + 5*v5 + 6*v6 + 7*v7 + 8*v8

    print(fit_multi_tex(eupper=u.Quantity(nacltbl['EU_K'], u.K),
                        nupperoverg=nacl_nu,
                        vstate=vstate,
                        jstate=np.array(nacltbl['J$_u$'], dtype='int'),
                        vibenergies=get_vib_energies(salt_tables.NaCl),
                        rotenergies=get_rot_energies(salt_tables.NaCl),
                        errors=enacl_nu,
                        rottemlims=(20,175),
                        vibtemlims=(500,8000),
                        collims=(np.log(1e10), np.log(1e16)),
                        plot=True, verbose=True, partition_func=partfunc, marker='o',
                        molname='NaCl',
                        colors=('r','g','b','orange','m','c','darkred','darkgreen',
                                'purple'), )
         )
    pl.legend(loc='upper right')
    pl.axis([300,4000,8.8,10.8])
    pl.title("NaCl")

    pl.savefig(paths.fpath("NaCl_rotational-vibrational_fit_diagrams.pdf"))







    #nacl37 = Vamdc.query_molecule('NaCl-37')
    #rt_nacl37 = nacl37.data['RadiativeTransitions']
    #frqs = u.Quantity([(float(rt_nacl37[key].FrequencyValue)*u.MHz).to(u.GHz,
    #                                                                   u.spectral())
    #                   for key in rt_nacl37])



    #upperStateRefs = [rt_nacl37[key].UpperStateRef for key in rt_nacl37]
    #degeneracies = [int(nacl37.data['States'][upperStateRef].TotalStatisticalWeight)
    #                for upperStateRef in upperStateRefs]
    #einsteinAij = u.Quantity([float(rt_nacl37[key].TransitionProbabilityA) for key in rt_nacl37], 1/u.s)

    frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('NaCl-37, v=0-15', catalog='CDMS', fmin=85*u.GHz, fmax=360*u.GHz)

    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    nacl37mask = np.array([(not hasattr(row['Species'], 'mask')) and
                          ('Na37Cl' == row['Species'][:6] or
                           '23Na-37Cl' in row['Species'])
                          for row in tbl])

    bad = (tbl['Line Name'] == 'Na37Clv=6_8-7') # tabulated velocity too far off

    print("Na37Cl: {0} in-band, {1} detected".format(nacl37mask.sum(),
                                                     (nacl37mask & (~bad)).sum()))

    nacl37tbl = tbl[nacl37mask & (~bad)]
    nacl37freqs = u.Quantity(nacl37tbl['Frequency'], u.GHz)
    kkms_nacl37 = (2*np.pi*nacl37tbl['Fitted Width']**2)**0.5 * nacl37tbl['Fitted Amplitude K']
    ekkms_nacl37 = (2*np.pi)**0.5 * (nacl37tbl['Fitted Width error']**2 *
                                     nacl37tbl['Fitted Amplitude K']**2 +
                                     nacl37tbl['Fitted Width']**2 *
                                     nacl37tbl['Fitted Amplitude error K']**2)**0.5
    #Aul = u.Quantity(list(map(get_Aul_(frqs), nacl37freqs)), u.Hz)
    #deg = u.Quantity(list(map(get_deg_(frqs), nacl37freqs)))
    Aul = nacl37tbl['Aij']
    deg = nacl37tbl['deg']
    nacl37_nu = nupper_of_kkms(kkms=kkms_nacl37,
                               freq=nacl37freqs,
                               Aul=Aul)/deg
    enacl37_nu = nupper_of_kkms(kkms=ekkms_nacl37,
                                freq=nacl37freqs,
                                Aul=Aul)/deg


    v0 = np.array(['v=0' in row['Species'] for row in nacl37tbl])
    v1 = np.array(['v=1' in row['Species'] for row in nacl37tbl])
    v2 = np.array(['v=2' in row['Species'] for row in nacl37tbl])
    v3 = np.array(['v=3' in row['Species'] for row in nacl37tbl])
    v4 = np.array(['v=4' in row['Species'] for row in nacl37tbl])
    v5 = np.array(['v=5' in row['Species'] for row in nacl37tbl])
    v6 = np.array(['v=6' in row['Species'] for row in nacl37tbl])
    j7 = np.array(['J=7' in row['Species'] for row in nacl37tbl])
    j8 = np.array(['J=8' in row['Species'] for row in nacl37tbl])

    pl.figure(4).clf()
    print("Na37Cl")
    tex0 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v0], u.K), nacl37_nu[v0],
                   errors=enacl37_nu[v0], plot=True,
                   verbose=True, partition_func=partfunc, marker='o', color='r')
    tex1 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v1], u.K), nacl37_nu[v1],
                   errors=enacl37_nu[v1], plot=True,
                   verbose=True, partition_func=partfunc, marker='s', color='b', label='v=1 ')
    tex2 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v2], u.K), nacl37_nu[v2], plot=True,
                   verbose=True, partition_func=partfunc, marker='^', color='g', label='v=2 ')
    #tex3 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v3], u.K), nacl37_nu[v3], plot=True,
    #               verbose=True, partition_func=partfunc, marker='o', color='r', label='v=3 ')
    #tex4 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v4], u.K), nacl37_nu[v4], plot=True,
    #               verbose=True, partition_func=partfunc, marker='d', color='orange', label='v=4 ')
    #tex5 = fit_tex(u.Quantity(nacl37tbl['EU_K'][v5], u.K), nacl37_nu[v5],
    #               errors=enacl37_nu[v5], plot=True,
    #               verbose=True, partition_func=partfunc, marker='v', color='m', label='v=5 ')
    texJ8 = fit_tex(u.Quantity(nacl37tbl['EU_K'][j8], u.K), nacl37_nu[j8],
                    errors=enacl37_nu[j8], plot=True, verbose=True,
                    partition_func=partfunc, marker='h', color='k', label='J=8 ')
    texJ7 = fit_tex(u.Quantity(nacl37tbl['EU_K'][j7], u.K), nacl37_nu[j7],
                    errors=enacl37_nu[j7], plot=True, verbose=True,
                    partition_func=partfunc, marker='>', color='c', label='J=7 ')

    pl.legend(loc='lower right')
    pl.title("Na$^{37}$Cl")
    pl.axis([-50,3000,9.0,12.25])
    pl.savefig(paths.fpath("NaCl37_rotational_diagrams.pdf"))


    pl.figure(4).clf()
    vstate = 0*v0 + 1*v1 + 2*v2 + 3*v3 + 4*v4 + 5*v5 + 6*v6

    print(fit_multi_tex(eupper=u.Quantity(nacl37tbl['EU_K'], u.K),
                        nupperoverg=nacl37_nu,
                        vstate=vstate,
                        jstate=np.array(nacl37tbl['J$_u$'], dtype='int'),
                        vibenergies=get_vib_energies(salt_tables.Na37Cl),
                        rotenergies=get_rot_energies(salt_tables.Na37Cl),
                        errors=enacl37_nu,
                        rottemlims=(20,225),
                        plot=True, verbose=True, partition_func=partfunc, marker='o',
                        molname='Na37Cl',
                        colors=('r','g','b','orange','m','c'), )
         )
    pl.legend(loc='lower right')
    pl.title("Na$^{37}$Cl")
    pl.axis((-50,3700,8,11.0))

    pl.savefig(paths.fpath("NaCl37_rotational-vibrational_fit_diagrams.pdf"))




    #k41cl = Vamdc.query_molecule('K-41-Cl')
    #rt41 = k41cl.data['RadiativeTransitions']
    #frqs = u.Quantity([(float(rt41[key].FrequencyValue)*u.MHz).to(u.GHz,
    #                                                              u.spectral())
    #                   for key in rt41])



    #upperStateRefs = [rt41[key].UpperStateRef for key in rt41]
    #degeneracies = [int(k41cl.data['States'][upperStateRef].TotalStatisticalWeight)
    #                for upperStateRef in upperStateRefs]
    #einsteinAij = u.Quantity([float(rt41[key].TransitionProbabilityA) for key in rt41], 1/u.s)

    frqs, einsteinAij, degeneracies, EU, partfunc = get_molecular_parameters('K-41-Cl, v=0-15', catalog='CDMS', fmin=85*u.GHz, fmax=360*u.GHz)

    tbl = table.Table.read(paths.tpath('fitted_stacked_lines.txt'), format='ascii.fixed_width')

    k41clmask = np.array([(not hasattr(row['Species'], 'mask')) and
                          ('41KCl' == row['Species'][:5] or
                           '41K-35Cl' in row['Species'])
                          for row in tbl])

    bad = (((tbl['Species'] == '41KClv=2') & (tbl['QNs'] == '29-28')) | # absorption
           ((tbl['Species'] == '41KClv=0') & (tbl['QNs'] == '31-30')) | # wing of NaCl
           ((tbl['Species'] == '41KClv=0') & (tbl['QNs'] == '46-45')) | # absorption
           (tbl['Line Name'] == 'NaClv=4_7-6') |
           (tbl['Line Name'] == '41KClv=2_29-28') |
           (tbl['Line Name'] == '41KClv=2 29-28') |
           (tbl['Line Name'] == '41KClv=0_31-30') |
           (tbl['Line Name'] == '41KClv=0 31-30') |
           (tbl['Line Name'] == '41KClv=0_46-45')
          )


    print("41KCl: {0} in-band, {1} detected".format(k41clmask.sum(),
                                                    (k41clmask & (~bad)).sum()))

    k41cltbl = tbl[k41clmask & (~bad)]
    k41clfreqs = u.Quantity(k41cltbl['Frequency'], u.GHz)
    kkms_k41cl = (2*np.pi*k41cltbl['Fitted Width']**2)**0.5 * k41cltbl['Fitted Amplitude K']
    ekkms_k41cl = (2*np.pi)**0.5 * (k41cltbl['Fitted Width error']**2 *
                                    k41cltbl['Fitted Amplitude K']**2 +
                                    k41cltbl['Fitted Width']**2 *
                                    k41cltbl['Fitted Amplitude error K']**2)**0.5
    #Aul = u.Quantity(list(map(get_Aul_(frqs), k41clfreqs)), u.Hz)
    #deg = u.Quantity(list(map(get_deg_(frqs), k41clfreqs)))
    Aul = k41cltbl['Aij']
    deg = k41cltbl['deg']
    k41cl_nu = nupper_of_kkms(kkms=kkms_k41cl,
                              freq=k41clfreqs,
                              Aul=Aul)/deg
    ek41cl_nu = nupper_of_kkms(kkms=ekkms_k41cl,
                               freq=k41clfreqs,
                               Aul=Aul)/deg

    v0 = np.array(['v=0' in row['Species'] for row in k41cltbl])
    v1 = np.array(['v=1' in row['Species'] for row in k41cltbl])
    v2 = np.array(['v=2' in row['Species'] for row in k41cltbl])
    v3 = np.array(['v=3' in row['Species'] for row in k41cltbl])

    pl.figure(5).clf()
    print("41KCl")
    tex0 = fit_tex(u.Quantity(k41cltbl['EU_K'][v0], u.K), k41cl_nu[v0],
                   errors=ek41cl_nu[v0], plot=True,
                   verbose=True, partition_func=partfunc, marker='o', color='r', label='v=0 ')
    #tex1 = fit_tex(u.Quantity(k41cltbl['EU_K'][v1], u.K), k41cl_nu[v1],
    #               errors=ek41cl_nu[v1], plot=True,
    #               verbose=True, partition_func=partfunc, marker='s', color='b', label='v=1 ')
    pl.legend(loc='lower right')
    pl.title("$^{41}$KCl")
    pl.savefig(paths.fpath("K41Cl_rotational_diagrams.pdf"))

    pl.figure(5).clf()
    vstate = 0*v0 + 1*v1 + 2*v2 + 3*v3
    print(fit_multi_tex(eupper=u.Quantity(k41cltbl['EU_K'], u.K),
                        nupperoverg=k41cl_nu,
                        vstate=vstate,
                        jstate=np.array(k41cltbl['J$_u$'], dtype='int'),
                        vibenergies=get_vib_energies(salt_tables.K41Cl),
                        rotenergies=get_rot_energies(salt_tables.K41Cl),
                        errors=ek41cl_nu,
                        plot=True, verbose=True, partition_func=partfunc, marker='o',
                        molname='41KCl',
                        colors=('r','g','b','orange','m','c',), )
         )
    pl.legend(loc='lower right')
    pl.title("$^{41}$KCl")

    pl.savefig(paths.fpath("41KCl_rotational-vibrational_fit_diagrams.pdf"))

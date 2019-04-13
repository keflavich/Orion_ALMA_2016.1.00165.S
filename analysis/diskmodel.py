# attempt to model Source I disk
import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from astropy import convolution
from astropy.modeling.functional_models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
import radio_beam
import lmfit
from constants import d_orion
from image_registration.fft_tools import shift
import paths
import regions
from astropy.nddata import Cutout2D
from mpl_plot_templates import asinh_norm
import latex_info
from latex_info import strip_trailing_zeros, round_to_n
from files import b7_hires_cont, b6_hires_cont, b3_hires_cont, reid7mm

import pylab as pl

# fh = fits.open('/Users/adam/work/orion/alma_lb/FITS/uid_A001_X88e_X1d3_calibrated_final_cont.pbcor.fits')
#
# cropslice_x = slice(3510,3610)
# cropslice_y = slice(3520,3620)
# data = fh[0].data.squeeze()[cropslice_y, cropslice_x]
# mywcs = wcs.WCS(fh[0].header).celestial[cropslice_y, cropslice_x]
# hdu = fits.PrimaryHDU(data=data, header=mywcs.to_header())
# hdu.writeto(paths.dpath('OrionSourceI_continuum_cutout_for_modeling.fits'),
#             overwrite=True)

epsfcn = 1e-3

fit_results = {}
ptsrc_diskcen_offsets = {}

parhist = {}

beams = {}

for fn, freq, band, thresh in [#('Orion_SourceI_B6_continuum_r-2_longbaselines_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                               #('Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                               #('Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                               (b7_hires_cont, 340.0*u.GHz, 'B7', 1*u.mJy),
                               (b6_hires_cont, 224.0*u.GHz, 'B6', 2*u.mJy),
                               (b3_hires_cont, 93.3*u.GHz, 'B3', None),
                               (reid7mm, 43.165*u.GHz, '7mm', None),
                              ]:
    print(band, freq)
    parhist[band] = {}

    fh = fits.open(paths.dpath(fn))

    coord = coordinates.SkyCoord("5:35:14.519", "-5:22:30.633", frame='fk5',
                                 unit=(u.hour, u.deg))
    mywcs = wcs.WCS(fh[0].header).celestial
    cutout = Cutout2D(data=fh[0].data.squeeze(),
                      wcs=mywcs,
                      position=coord,
                      size=0.5*u.arcsec)

    data = cutout.data
    new_header = cutout.wcs.to_header()
    mywcs = cutout.wcs

    fits.PrimaryHDU(data=data,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_fitted_data.fits'
                            .format(band)),
                        overwrite=True)

    pixscale = wcs.utils.proj_plane_pixel_area(mywcs)**0.5 * u.deg

    # old version diskends = coordinates.SkyCoord(['5:35:14.5232 -5:22:30.73',
    # old version                                  '5:35:14.5132 -5:22:30.54'],
    # old version                                 frame='fk5', unit=(u.hour, u.deg))
    diskend_regs = regions.read_ds9(paths.rpath('diskends.reg'))
    diskends = coordinates.SkyCoord([reg.center for reg in diskend_regs])

    diskends_pix = np.array(mywcs.wcs_world2pix(diskends.ra.deg, diskends.dec.deg, 0))
    (x1,x2),(y1,y2) = diskends_pix

    observed_beam = radio_beam.Beam.from_fits_header(fh[0].header)
    beams[band] = observed_beam

    data_K = (data*u.Jy).to(u.K, observed_beam.jtok_equiv(freq))

    new_header_K = new_header.copy()
    new_header_K['BUNIT'] = 'K'
    fits.PrimaryHDU(data=data_K.value, header=new_header_K).writeto(
        paths.dpath('contmodels/{0}_data_K_cutout.fits'.format(band)),
        overwrite=True)

    jtok = observed_beam.jtok(freq)

    # from the first round of fitting, there is a residual source at this position
    ptsrc = coordinates.SkyCoord(83.81049240934931*u.deg, -5.375170355557261*u.deg, frame='icrs')

    ptsrcx, ptsrcy = mywcs.wcs_world2pix(ptsrc.ra, ptsrc.dec, 0)
    ptsrc_amp_value = data[int(ptsrcy), int(ptsrcx)]
    print("Point source amplitude data value: {0}".format(ptsrc_amp_value))
    print()

    if thresh is not None:
        # mask out low pixels: see what happens if we only fit the stuff we
        # *really* believe.
        baddata = u.Quantity(data, u.Jy) < thresh
        ndata = (~baddata).sum()
    else:
        baddata = None
        ndata = data.size

    print("Threshold: {0} ndata: {1}".format(thresh, ndata))

    def model(x1, x2, y1, y2, scale, kernelmajor=None, kernelminor=None, kernelpa=None,
              ptsrcx=None, ptsrcy=None, ptsrcamp=None, ptsrcwid=None):
        """
        The model, with a variable number of parameters....
        """
        #if any(ii < 0 for ii in (x1,x2,y1,y2)):
        #    return 1e5
        #if (x1 > data.shape[1]-1) or (x2 > data.shape[1]-1) or (y1 > data.shape[0]-1) or (y2 > data.shape[0]-1):
        #    return 1e5
        if x2 > data.shape[1] - 1:
            x2 = data.shape[1] - 1
        if y2 > data.shape[0] - 1:
            y2 = data.shape[0] - 1

        def line_y(x):
            # y = m x + b
            m = (y2-y1)/(x2-x1)
            b = y1 - x1*m
            return m*x+b

        xx = np.linspace(x1, x2, 1000)
        yy = line_y(xx)

        if hasattr(scale, 'value'):
            scale = scale.value

        disk = np.zeros_like(data, dtype='float')
        disk[(np.round(yy).astype('int')), (np.round(xx).astype('int'))] = scale

        if kernelmajor is None:
            beam = observed_beam
        else:
            beam = radio_beam.Beam(kernelmajor*u.arcsec, kernelminor*u.arcsec,
                                   kernelpa*u.deg).convolve(observed_beam)

        beam_kernel = beam.as_kernel(pixscale).array
        # want a peak-normalized kernel
        beam_amp = beam_kernel.max()

        diskmod = convolution.convolve_fft(disk, beam_kernel) / beam_amp

        if ptsrcwid is not None:
            assert ptsrcamp is not None

            #posang = np.arctan2(y2 - y1, x2 - x1)*u.rad - 90*u.deg

            ptsrc_wid_bm = radio_beam.Beam(ptsrcwid*u.arcsec, 0.00001*u.arcsec,
                                           (kernelpa+90)*u.deg)
            convbm = observed_beam.convolve(ptsrc_wid_bm)
            assert convbm.pa.value != 0
            ptsrc_bm = convbm.as_kernel(pixscale,
                                        x_size=data.shape[1],
                                        y_size=data.shape[0])
            ptsrcmod = (shift.shift2d(ptsrc_bm,
                                      ptsrcx-data.shape[0]/2, ptsrcy-data.shape[1]/2) /
                        beam_amp * ptsrcamp)
            diskmod += ptsrcmod
        elif ptsrcamp is not None:
            ptsrcmod = (shift.shift2d(observed_beam.as_kernel(pixscale,
                                                              x_size=data.shape[1],
                                                              y_size=data.shape[0]),
                                      ptsrcx-data.shape[0]/2, ptsrcy-data.shape[1]/2) /
                        beam_amp * ptsrcamp)
            diskmod += ptsrcmod


        return diskmod

    obsbm_max = observed_beam.as_kernel(pixscale).array.max()

    def ptsrcmodel(ptsrcx, ptsrcy, ptsrcamp, ptsrcwid, kernelpa=52):
        ptsrc_wid_bm = radio_beam.Beam(ptsrcwid*u.arcsec, 0.00001*u.arcsec,
                                       (kernelpa+90)*u.deg)
        convbm = observed_beam.convolve(ptsrc_wid_bm)
        assert convbm.pa.value != 0
        ptsrc_bm = convbm.as_kernel(pixscale,
                                    x_size=data.shape[1],
                                    y_size=data.shape[0])
        ptsrcmod = (shift.shift2d(ptsrc_bm,
                                  ptsrcx-data.shape[0]/2, ptsrcy-data.shape[1]/2) /
                    obsbm_max * ptsrcamp)
        return ptsrcmod

    def residual(pars, rms=0.0001, maskptsrc=False, data=data, model=model):
        mod = model(**pars)

        if maskptsrc:
            yy,xx = np.indices(data.shape)
            rr = ((yy-ptsrcy)**2 + (xx-ptsrcx)**2)**0.5
            mask = rr > 5. # 5 is arbitrary...
            if baddata is not None:
                mask &= ~baddata
        else:
            if baddata is not None:
                mask = ~baddata
            else:
                mask = None

        if mask is None:
            # added abs val just to make extra sure...
            resid = ((data - mod)/rms)
        else:
            resid = ((data[mask] - mod[mask])/rms)

        parhist[band][len(pars)].append([pars.valuesdict(), (resid**2).sum()])

        assert (resid**2).sum() == (np.abs(resid)**2).sum()
        return resid

    diskmod = model(x1,x2,y1,y2,data.max())
    print("Initial disk model max: {0} data max: {1}".format(diskmod.max(), data.max()))

    parameters = lmfit.Parameters()
    parameters.add('x1', value=x1)
    parameters.add('x2', value=x2)
    parameters.add('y1', value=y1)
    parameters.add('y2', value=y2)
    parameters.add('scale', value=data.max())
    parhist[band][len(parameters)] = []
    result = lmfit.minimize(residual, parameters, epsfcn=epsfcn, kws={'maskptsrc':False})
    print("Basic fit parameters (linear model):")
    result.params.pretty_print()
    #print("red Chi^2: {0:0.3g}".format(result.chisqr / (ndata - result.nvarys)))
    print("red Chi^2: {0:0.3g}".format(result.redchi))
    print(result.message)
    print()

    bestdiskmod_beam = model(**result.params)

    # Create a "beam" that is really the vertical x horizontal scale height
    # to be convolved with the observed beam
    # old versions for QA2 cutout parameters.add('kernelmajor', value=0.064)
    # old versions for QA2 cutout parameters.add('kernelminor', value=0.043)
    parameters.add('kernelmajor', value=0.054, min=0.01, max=0.3)
    parameters.add('kernelminor', value=0.033, min=0.01, max=0.3)
    # Fix the position angle such that one direction of the resulting kernel will
    # directly be a Gaussian scale height
    measured_positionangle = -38
    parameters.add('kernelpa', value=measured_positionangle+90, vary=False)
    parhist[band][len(parameters)] = []
    result2 = lmfit.minimize(residual, parameters, epsfcn=epsfcn, kws={'maskptsrc':False})
    print("Smoothed linear fit parameters:")
    result2.params.pretty_print()
    #print("red Chi^2: {0:0.3g}".format(result2.chisqr / (ndata - result2.nvarys)))
    print("red Chi^2: {0:0.3g}".format(result2.redchi))
    print(result2.message)
    print()

    bestdiskmod = model(**result2.params)

    parameters.add('ptsrcx', value=ptsrcx, min=ptsrcx-5, max=ptsrcx+5)
    parameters.add('ptsrcy', value=ptsrcy, min=ptsrcy-5, max=ptsrcy+5)
    parameters.add('ptsrcamp', value=0.004, min=0.001, max=0.6)
    parhist[band][len(parameters)] = []
    result3 = lmfit.minimize(residual, parameters, epsfcn=epsfcn)
    print("Smoothed linear fit parameters with point source:")
    result3.params.pretty_print()
    #print("red Chi^2: {0:0.3g}".format(result3.chisqr / (ndata - result3.nvarys)))
    print("red Chi^2: {0:0.3g}".format(result3.redchi))
    print(result3.message)
    print()

    bestdiskplussourcemod = model(**result3.params)

    parameters.add('ptsrcwid', value=0.04, min=0.01, max=0.1)
    parhist[band][len(parameters)] = []
    minimizer = lmfit.Minimizer(residual, parameters, epsfcn=epsfcn)
    result4 = minimizer.minimize()
    print("Smoothed linear fit parameters with horizontally smeared point source:")
    result4.params.pretty_print()
    #print("red Chi^2: {0:0.3g}".format(result4.chisqr / (ndata - result4.nvarys)))
    print("red Chi^2: {0:0.3g}".format(result4.redchi))
    print(result4.message)
    print()

    bestdiskplussmearedsourcemod = model(**result4.params)


    # remove the source
    sourceless_pars = dict(**result4.params)
    sourceless_pars['ptsrcamp'] = 0
    bestdiskminussmearedsourcemod = model(**sourceless_pars)

    # refit with source position, amplitude as only free parameters in the
    # hope that we can estimate errors
    source_pars = result4.params.copy()
    for pn in ('x1','x2','y1','y2','scale', 'kernelmajor', 'kernelminor', 'kernelpa'):
        source_pars.pop(pn)
    source_pars['ptsrcwid'].vary=False

    data_nodisk = data - bestdiskminussmearedsourcemod
    data_noptsrc = data - (bestdiskplussmearedsourcemod-bestdiskminussmearedsourcemod)

    parhist[band][len(source_pars)] = []
    result4b = lmfit.minimize(residual, source_pars, epsfcn=epsfcn,
                              kws=dict(data=data_nodisk, model=ptsrcmodel))
    print("Smoothed linear fit parameters with horizontally smeared point source (many pars held fixed):")
    result4b.params.pretty_print()
    #print("red Chi^2: {0:0.3g}".format(result4.chisqr / (ndata - result4.nvarys)))
    print("red Chi^2: {0:0.3g}".format(result4b.redchi))
    print(result4b.message)
    print()

    y_, x_ = np.indices(data_nodisk.shape)
    apylmfit = LevMarLSQFitter()
    # the weights don't matter.
    # I feel like they should, but maybe there's some theorem I'm not
    # appreciating that says that the parameter constraints are independent of
    # the errors on the data...
    # if I'm wrong, this is a bug in astropy.
    astropy_fit_results = apylmfit(Gaussian2D(amplitude=source_pars['ptsrcamp'].value,
                                              x_mean=source_pars['ptsrcx'].value,
                                              y_mean=source_pars['ptsrcy'].value,
                                              x_stddev=2,),
                                   x_, y_, data_nodisk,
                                   weights=1/(np.ones_like(data_nodisk) *
                                              (data-bestdiskminussmearedsourcemod).std()),
                                  )

    print("astropy gaussian fit to residual point source image: ")
    print(astropy_fit_results)
    print(np.diagonal(apylmfit.fit_info['param_cov'])**0.5)
    astropy_positional_offset = (apylmfit.fit_info['param_cov'][1,1] +
                                 apylmfit.fit_info['param_cov'][2,2])**0.5
    print(astropy_positional_offset)
    print()



    refit_sourceless_params = result4.params.copy()
    refit_sourceless_params.pop('ptsrcx')
    refit_sourceless_params.pop('ptsrcy')
    refit_sourceless_params.pop('ptsrcamp')
    refit_sourceless_params.pop('ptsrcwid')

    parhist[band][len(refit_sourceless_params)] = []
    result6 = lmfit.minimize(residual, refit_sourceless_params, epsfcn=epsfcn)
    print("Smoothed linear fit parameters without horizontally smeared point source, refit:")
    result6.params.pretty_print()
    print("red Chi^2: {0:0.3g}".format(result6.redchi))
    print(result6.message)
    print()

    bestdiskmod_refit = model(**result6.params)

    print("For comparison, the original smoothed linear fit: ")
    result2.params.pretty_print()
    print("red Chi^2: {0:0.3g}".format(result2.redchi))
    print()

    parameters.pop('kernelmajor')
    parameters.pop('kernelminor')
    parhist[band][len(parameters)] = []
    result5 = lmfit.minimize(residual, parameters, epsfcn=epsfcn)
    print("Linear fit parameters with horizontally smeared point source:")
    result5.params.pretty_print()
    print("red Chi^2: {0:0.3g}".format(result5.redchi))
    print(result5.message)
    print()

    bestbasicdiskplussmearedsourcemod = model(**result5.params)


    parhist[band][len(refit_sourceless_params)] = []
    result7 = lmfit.minimize(residual, refit_sourceless_params, epsfcn=epsfcn,
                             kws={'data':data_noptsrc})
    print("Smoothed linear fit parameters with best point source removed from data:")
    result7.params.pretty_print()
    print("red Chi^2: {0:0.3g}".format(result7.redchi))
    print(result7.message)
    print()

    bestdiskmod_refit_withoutptsrc = model(**result7.params)


    fits.PrimaryHDU(data=bestbasicdiskplussmearedsourcemod,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_bestbasicdiskplussmearedsourcemodel.fits'
                            .format(band)),
                        overwrite=True)
    fits.PrimaryHDU(data=bestdiskplussmearedsourcemod,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_bestdiskplussmearedsourcemodel.fits'
                            .format(band)),
                        overwrite=True)
    fits.PrimaryHDU(data=bestdiskplussourcemod,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_bestdiskplussourcemodel.fits'
                            .format(band)),
                        overwrite=True)
    fits.PrimaryHDU(data=bestdiskmod,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_bestdiskmodel.fits'
                            .format(band)),
                        overwrite=True)
    fits.PrimaryHDU(data=bestdiskmod_beam,
                    header=new_header).writeto(
                        paths.dpath(
                            'contmodels/{0}_bestdiskmod_beamel.fits'
                            .format(band)),
                        overwrite=True)


    ptsrc_ra, ptsrc_dec = mywcs.wcs_pix2world(result4.params['ptsrcx'],
                                              result4.params['ptsrcy'], 0)
    try:
        eptsrc_ra, eptsrc_dec = result4.params['ptsrcx'].stderr*pixscale, result4.params['ptsrcy'].stderr*pixscale,
    except TypeError:
        eptsrc_ra, eptsrc_dec = np.nan*u.arcsec,np.nan*u.arcsec
    fitted_ptsrc = coordinates.SkyCoord(ptsrc_ra*u.deg, ptsrc_dec*u.deg, frame=mywcs.wcs.radesys.lower())
    print("Fitted point source location = {0} {1}".format(fitted_ptsrc.to_string('hmsdms'), fitted_ptsrc.frame.name))
    print("Fitted point source amplitude: {0}".format(result4.params['ptsrcamp']))

    print("diskends: {0}".format(diskends))
    fitted_diskends_mod1 = coordinates.SkyCoord(*mywcs.wcs_pix2world([result.params['x1'], result.params['x2']], [result.params['y1'], result.params['y2']], 0),
                                                unit=(u.deg, u.deg),
                                                frame=mywcs.wcs.radesys.lower())
    print("fitted diskends (model 1): {0}".format(fitted_diskends_mod1))
    fitted_diskends_mod2 = coordinates.SkyCoord(*mywcs.wcs_pix2world([result2.params['x1'], result2.params['x2']], [result2.params['y1'], result2.params['y2']], 0),
                                                unit=(u.deg, u.deg),
                                                frame=mywcs.wcs.radesys.lower())
    print("fitted diskends (model 2): {0}".format(fitted_diskends_mod2))
    fitted_diskends_mod3 = coordinates.SkyCoord(*mywcs.wcs_pix2world([result3.params['x1'], result3.params['x2']], [result3.params['y1'], result3.params['y2']], 0),
                                                unit=(u.deg, u.deg),
                                                frame=mywcs.wcs.radesys.lower())
    print("fitted diskends (model 3): {0}".format(fitted_diskends_mod3))
    fitted_diskends_mod4 = coordinates.SkyCoord(*mywcs.wcs_pix2world([result4.params['x1'], result4.params['x2']], [result4.params['y1'], result4.params['y2']], 0),
                                                unit=(u.deg, u.deg),
                                                frame=mywcs.wcs.radesys.lower())
    print("fitted diskends (model 4): {0}".format(fitted_diskends_mod4))

    disk_center_mod4 = coordinates.SkyCoord(fitted_diskends_mod4.ra.mean(),
                                            fitted_diskends_mod4.dec.mean(),
                                            frame=fitted_diskends_mod4.frame.name)

    ptsrc_diskcen_sep = fitted_ptsrc.separation(disk_center_mod4).to(u.arcsec)
    print("fitted pointsource is offset from center by {0}".format(ptsrc_diskcen_sep))
    ptsrc_diskcen_offsets[band] = ptsrc_diskcen_sep, astropy_positional_offset, astropy_positional_offset*pixscale

    posang = np.arctan2(result3.params['y2']-result3.params['y1'],
                        result3.params['x2']-result3.params['x1'])*u.rad - 90*u.deg
    print("posang={0}".format(90*u.deg+posang.to(u.deg)))
    print("kernelpa2={0}".format(result2.params['kernelpa']))
    print("kernelpa3={0}".format(result3.params['kernelpa']))
    print("kernelpa4={0}".format(result4.params['kernelpa']))
    print("kernelpa5={0}".format(result5.params['kernelpa']))

    print()
    print("Disk center is {0} in model 4".format(disk_center_mod4))
    ll_0p4 = coordinates.SkyCoord(disk_center_mod4.ra - 0.4*u.arcsec/np.cos(disk_center_mod4.dec) * np.sin(posang),
                                  disk_center_mod4.dec - 0.4*u.arcsec * np.cos(posang),
                                  frame=fitted_diskends_mod4.frame.name)
    ur_0p4 = coordinates.SkyCoord(disk_center_mod4.ra + 0.4*u.arcsec/np.cos(disk_center_mod4.dec) * np.sin(posang),
                                  disk_center_mod4.dec + 0.4*u.arcsec * np.cos(posang),
                                  frame=fitted_diskends_mod4.frame.name)
    print("lower left, upper right corners are: {0}, {1}".format(ll_0p4.icrs, ur_0p4.icrs))
    print("line({0},{1},{2},{3})".format(ll_0p4.icrs.ra.deg,
                                         ll_0p4.icrs.dec.deg,
                                         ur_0p4.icrs.ra.deg,
                                         ur_0p4.icrs.dec.deg))
    print("line({0},{1},{2},{3})".format(fitted_diskends_mod4.icrs.ra.deg[0],
                                         fitted_diskends_mod4.icrs.dec.deg[0],
                                         fitted_diskends_mod4.icrs.ra.deg[1],
                                         fitted_diskends_mod4.icrs.dec.deg[1],))

    # compute offset from point source to disk center along the disk axis angle
    assert ptsrc_diskcen_sep < 0.1*u.arcsec
    sep_projected = ptsrc_diskcen_sep * np.abs(np.cos(posang-90*u.deg))
    print("The source is separated from the disk center by {0} = {1} in projection"
          .format(sep_projected, (sep_projected*d_orion).to(u.au,
                                                            u.dimensionless_angles())))

    pointmodel_image = model(0,1,0,1,scale=0,
                             kernelmajor=result4.params['kernelmajor'],
                             kernelminor=result4.params['kernelminor'],
                             kernelpa=result4.params['kernelpa'],
                             ptsrcx=result4.params['ptsrcx'],
                             ptsrcy=result4.params['ptsrcy'],
                             ptsrcamp=result4.params['ptsrcamp'],
                             ptsrcwid=result4.params['ptsrcwid'],
                            )

    fitted_beam = radio_beam.Beam(result2.params['kernelmajor']*u.arcsec,
                                  result2.params['kernelminor']*u.arcsec,
                                  result2.params['kernelpa']*u.deg,)
    # NOTE: the fitted beam *is* the source size after a revision to the model
    # in which the input beam is convolved with the observed beam
    #source_size = fitted_beam.deconvolve(observed_beam)
    source_size = fitted_beam
    print("Result2 version (no point source): ")
    print("Fitted source (disk vertical scale) size: {0}".format(fitted_beam.__repr__()))
    print("Real source (disk vertical scale) size: {0}".format(source_size.__repr__()))
    scaleheight = (fitted_beam.major*d_orion).to(u.au, u.dimensionless_angles())
    print("Scale height: {0}".format(scaleheight))
    print()

    # the fit is much better with the source included...
    fitted_beam = radio_beam.Beam(result4.params['kernelmajor']*u.arcsec,
                                  result4.params['kernelminor']*u.arcsec,
                                  result4.params['kernelpa']*u.deg,)
    source_size = fitted_beam
    print("Result4 version (has smeared point source): ")
    print("Fitted source (disk vertical scale) size: {0}".format(fitted_beam.__repr__()))
    print("Real source (disk vertical scale) size: {0}".format(source_size.__repr__()))
    scaleheight = (fitted_beam.major*d_orion).to(u.au, u.dimensionless_angles())
    try:
        scaleheight_error = (result4.params['kernelmajor'].stderr*u.arcsec*d_orion).to(u.au, u.dimensionless_angles())
    except TypeError:
        scaleheight_error = np.nan * u.au
    print("Scale height: {0} +/- {1}".format(scaleheight, scaleheight_error))
    print()


    length_as = (((result4.params['x2'] - result4.params['x1'])**2 +
                  (result4.params['y2'] - result4.params['y1'])**2)**0.5 * pixscale).to(u.arcsec)
    length_au = (length_as * d_orion).to(u.au, u.dimensionless_angles())

    # complicated propagation of errors using rules I didn't even know about...
    x1, ex1 = result.params['x1'].value, result.params['x1'].stderr
    x2, ex2 = result.params['x2'].value, result.params['x2'].stderr
    y1, ey1 = result.params['y1'].value, result.params['y1'].stderr
    y2, ey2 = result.params['y2'].value, result.params['y2'].stderr
    eA2 = 2*(x2**2*ex2**2 + ex2**2*x1**2 + ex1**2*x2**2 + x1**2*ex1**2)**0.5
    eB2 = 2*(y2**2*ey2**2 + ey2**2*y1**2 + ey1**2*y2**2 + y1**2*ey1**2)**0.5
    eW = (eA2**2+eB2**2)**0.5
    W = ((x2-x1)**2+(y2-y1)**2)
    elength = (0.5 * W**-0.5 * eW)
    elength_au = (elength * pixscale * d_orion).to(u.au, u.dimensionless_angles())

    print("Length in arcsec: {0:0.3g}  in AU: {1:0.3g} +/- {3:0.3g}  or radius {2:0.3g}"
          .format(length_as, length_au, length_au/2, elength_au))
    print()
    print()

    ppbeam = (observed_beam.sr/pixscale**2).decompose()

    print("pointmodel_image_sum/ppbeam = {0}  ptsrcamp = {1}".format(pointmodel_image.sum()/ppbeam,
                                                                     result4.params['ptsrcamp'].value))

    with open(paths.rpath('{0}_continuum_disk.reg'.format(band)), 'w') as fh:
        fh.write('icrs\n')
        fh.write("line({0},{1},{2},{3}) # text={{{4}}} color=magenta\n"
                 .format(fitted_diskends_mod4.icrs.ra.deg[0],
                         fitted_diskends_mod4.icrs.dec.deg[0],
                         fitted_diskends_mod4.icrs.ra.deg[1],
                         fitted_diskends_mod4.icrs.dec.deg[1],
                         band
                        ))
        fh.write("line({0},{1},{2},{3}) # text={{{4} pm0.4}} color=green\n"
                 .format(ll_0p4.icrs.ra.deg, ll_0p4.icrs.dec.deg,
                         ur_0p4.icrs.ra.deg, ur_0p4.icrs.dec.deg,
                         band
                        ))

        fh.write("point({0},{1}) # text={{{2} diskcen}} point=x color=red\n"
                 .format(disk_center_mod4.icrs.ra.deg,
                         disk_center_mod4.icrs.dec.deg,
                         band
                        ))
        fh.write("point({0},{1}) # text={{{2} ptsrc}} point=circle color=blue\n"
                 .format(fitted_ptsrc.icrs.ra.deg,
                         fitted_ptsrc.icrs.dec.deg,
                         band
                        ))

    try:
        ePtwidth = (u.Quantity(result4.params['ptsrcwid'].stderr, u.arcsec)*d_orion).to(u.au, u.dimensionless_angles()).value
    except TypeError:
        ePtwidth = np.nan*u.au

    fit_results[freq] = {
                         'Disk FWHM': scaleheight,
                         'eDisk FWHM': scaleheight_error,
                         'Disk Radius': length_au/2,
                         'eDisk Radius': elength_au,
                         'Disk PA': posang,
                         'Pt Position': fitted_ptsrc,
                         'ePt RA': eptsrc_ra.to(u.arcsec).value*15,
                         'ePt Dec': eptsrc_dec.to(u.arcsec).value,
                         'Pt Amp': result4.params['ptsrcamp'].value,
                         'ePt Amp': result4.params['ptsrcamp'].stderr,
                         'Pt Width': (u.Quantity(result4.params['ptsrcwid'], u.arcsec)*d_orion).to(u.au, u.dimensionless_angles()).value,
                         'ePt Width': ePtwidth,
                         'Pt Flux': pointmodel_image.sum() / ppbeam,
                         'Total Flux': bestdiskplussmearedsourcemod.sum() / ppbeam,
                         'Pt \%': pointmodel_image.sum() / bestdiskplussmearedsourcemod.sum(),
                        }


    pl.figure(1)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestdiskmod_beam, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestdiskmod_beam, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestdiskmod_beam, colors=['w']*100, levels=np.linspace(bestdiskmod_beam.max()*0.05, bestdiskmod_beam.max(), 5))
    axlims = pl.axis()
    pl.plot([result.params['x1'], result.params['x2']],
            [result.params['y1'], result.params['y2']],
            'r')
    pl.axis(axlims)

    pl.savefig(paths.fpath("contmodel/SourceI_Disk_model_{0}.pdf".format(band)), bbox_inches='tight')

    pl.figure(2)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestdiskmod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestdiskmod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestdiskmod, colors=['w']*100, levels=np.linspace(bestdiskmod.max()*0.05, bestdiskmod.max(), 5))
    axlims = pl.axis()
    pl.plot([result2.params['x1'], result2.params['x2']],
            [result2.params['y1'], result2.params['y2']],
            'r')
    pl.axis(axlims)

    pl.savefig(paths.fpath("contmodel/SourceI_Disk_model_bigbeam_{0}.pdf".format(band)), bbox_inches='tight')

    pl.figure(3)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestdiskplussourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestdiskplussourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestdiskplussourcemod, colors=['w']*100, levels=np.linspace(bestdiskplussourcemod.max()*0.05, bestdiskplussourcemod.max(), 5))
    axlims = pl.axis()
    pl.plot([result2.params['x1'], result2.params['x2']],
            [result2.params['y1'], result2.params['y2']],
            'r')
    pl.axis(axlims)

    pl.savefig(paths.fpath("contmodel/SourceI_Disk_model_bigbeam_withptsrc_{0}.pdf".format(band)), bbox_inches='tight')

    pl.figure(4)
    pl.clf()
    pl.imshow(data_K.value, interpolation='none', origin='lower', cmap='viridis')
    pl.colorbar()

    pl.figure(7)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestdiskplussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestdiskplussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestdiskplussmearedsourcemod, colors=['w']*100, levels=np.linspace(bestdiskplussmearedsourcemod.max()*0.05, bestdiskplussmearedsourcemod.max(), 5))

    ptsrc_wid_bm = radio_beam.Beam(result4.params['ptsrcwid']*u.arcsec, 0.00001*u.arcsec, observed_beam.pa-90*u.deg) #result4.params['kernelpa']*u.deg)
    ptsrc_wid_bm = radio_beam.Beam(result4.params['ptsrcwid']*u.arcsec, 0.00001*u.arcsec, 90*u.deg+result4.params['kernelpa']*u.deg)
    ptsrc_bm = observed_beam.convolve(ptsrc_wid_bm)
    for ax in pl.gcf().axes:
        if band == 'B6':
            ell_smearing_bm = ptsrc_wid_bm.ellipse_to_plot(20, 85, pixscale)
            ell_wid_bm = ptsrc_bm.ellipse_to_plot(20, 70, pixscale)
            ell_obs_bm = observed_beam.ellipse_to_plot(20, 100, pixscale)
        elif band == 'B3':
            ell_smearing_bm = ptsrc_wid_bm.ellipse_to_plot(10, 45, pixscale)
            ell_wid_bm = ptsrc_bm.ellipse_to_plot(10, 35, pixscale)
            ell_obs_bm = observed_beam.ellipse_to_plot(10, 55, pixscale)
        else:
            continue
        for ell in (ell_smearing_bm, ell_wid_bm, ell_obs_bm):
            ax.add_patch(ell)

    axlims = pl.axis()
    pl.plot([result2.params['x1'], result2.params['x2']],
            [result2.params['y1'], result2.params['y2']],
            'r')
    pl.axis(axlims)

    pl.savefig(paths.fpath("contmodel/SourceI_Disk_model_bigbeam_withsmearedptsrc_{0}.pdf".format(band)), bbox_inches='tight')


    pl.figure(8)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestbasicdiskplussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestbasicdiskplussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestbasicdiskplussmearedsourcemod, colors=['w']*100,
               levels=np.linspace(bestbasicdiskplussmearedsourcemod.max()*0.05,
                                  bestbasicdiskplussmearedsourcemod.max(), 5))
    pl.savefig(paths.fpath("contmodel/SourceI_BasicDisk_model_withsmearedptsrc_{0}.pdf".format(band)), bbox_inches='tight')


    pl.figure(9)
    pl.clf()
    pl.subplot(2,2,1)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    #pl.plot(diskends_pix[0,:], diskends_pix[1,:], 'w')
    #pl.plot(xx, yy, 'r')
    pl.subplot(2,2,2)
    pl.imshow(bestdiskminussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,3)
    pl.imshow(data - bestdiskminussmearedsourcemod, interpolation='none', origin='lower', cmap='viridis')
    pl.subplot(2,2,4)
    pl.imshow(data, interpolation='none', origin='lower', cmap='viridis')
    pl.contour(bestdiskminussmearedsourcemod, colors=['w']*100,
               levels=np.linspace(bestdiskminussmearedsourcemod.max()*0.05,
                                  bestdiskminussmearedsourcemod.max(), 5))
    pl.savefig(paths.fpath("contmodel/SourceI_SmoothDisk_model_withOUTsmearedptsrc_{0}.pdf".format(band)), bbox_inches='tight')





    # publication figures
    fig5 = pl.figure(5)
    fig5.clf()
    ax = fig5.gca()
    im0 = ax.imshow(data*jtok.value, interpolation='none', origin='lower', cmap='viridis')
    im = ax.imshow(data*1e3, interpolation='none', origin='lower', cmap='viridis')

    wavelength = '1.3 mm' if band == 'B6' else '3 mm' if band == 'B3' else '7 mm' if band == '7mm' else '880 \mu\mathrm{m}' if band == 'B7' else 'ERROR'

    cb = fig5.colorbar(mappable=im)
    cb.set_label('$S_{{{0}}}$ [mJy beam$^{{-1}}$]'.format(wavelength))
    cb.ax.set_aspect('auto')
    pos = cb.ax.get_position()
    cb2 = cb.ax.twinx()
    cb2.set_ylim([-5, 40])
    pos.x0 += 0.06
    cb.ax.set_position(pos)
    cb2.set_position(pos)
    cb.ax.yaxis.set_label_position("left")
    cb2.set_ylabel('$T_B$ [K]')


    ax.set_xticks([])
    ax.set_yticks([])

    beam_ellipse = observed_beam.ellipse_to_plot(0.9*data.shape[1], 0.1*data.shape[0], pixscale)
    beam_ellipse.set_edgecolor('w')
    beam_ellipse.set_facecolor('w')
    ax.add_patch(beam_ellipse)
    fig5.savefig(paths.fpath("contmodel/OrionSourceI_data_{0}.pdf".format(band)), bbox_inches='tight')

    extent = [-data.shape[1]/2*pixscale.to(u.arcsec).value,
              data.shape[1]/2*pixscale.to(u.arcsec).value,
              -data.shape[0]/2*pixscale.to(u.arcsec).value,
              data.shape[0]/2*pixscale.to(u.arcsec).value,]

    plotted_center = mywcs.wcs_pix2world(data.shape[1]/2, data.shape[0]/2, 0)
    plotted_center_coord = coordinates.SkyCoord(*plotted_center, unit=(u.deg,
                                                                       u.deg),
                                                frame=mywcs.wcs.radesys.lower())
    print("Center position of {2} image is: {0}  {1}"
          .format(plotted_center_coord.to_string('hmsdms'),
                  plotted_center_coord.to_string('hmsdms', sep=':'),
                  band
                 )
         )

    fig5.clf()
    ax = fig5.gca()
    #im0 = ax.imshow(data*jtok.value, interpolation='none', origin='lower',
    #                cmap='gray', vmin=-5, vmax=40,
    #                extent=extent,
    #               )
    im = ax.imshow(data*1e3, interpolation='none', origin='lower',
                   cmap='gray_r', vmin=-5/jtok.value*1e3,
                   vmax=40/jtok.value*1e3,
                   extent=extent,
                  )
    con = ax.contour(data*jtok.value,
                     levels=[50, 100, 150, 200, 300, 400, 500, 600],
                     colors=['w']*10,
                     linewidths=0.75,
                     extent=extent,
                    )
    #cb2 = fig5.colorbar(mappable=im0)
    # https://stackoverflow.com/questions/27151098/draw-colorbar-with-twin-scales

    ebm = observed_beam.ellipse_to_plot(extent[0]+0.05, extent[2]+0.05, 1./3600*u.deg)
    ebm.set_facecolor('none')
    ebm.set_edgecolor('b')
    ax.add_patch(ebm)

    cb = fig5.colorbar(mappable=im)
    cb.set_label('$S_{{{0}}}$ [mJy beam$^{{-1}}$]'.format(wavelength))
    cb.ax.set_aspect('auto')
    pos = cb.ax.get_position()
    cb2 = cb.ax.twinx()
    cb2.set_ylim([-5, 40])
    pos.x0 +=0.06
    cb.ax.set_position(pos)
    cb2.set_position(pos)
    cb.ax.yaxis.set_label_position("left")
    cb2.set_ylabel('$T_B$ [K]')


    #ax.set_xticks([])
    #ax.set_yticks([])
    ax.set_xlabel("RA offset (\")")
    ax.set_ylabel("Dec offset (\")")

    # need duplicate code here (grr) because matplotlib refuses to reuse Aritst objects
    beam_ellipse = observed_beam.ellipse_to_plot(0.9*data.shape[1], 0.1*data.shape[0], pixscale)
    beam_ellipse.set_edgecolor('w')
    beam_ellipse.set_facecolor('w')
    ax.add_patch(beam_ellipse)
    fig5.savefig(paths.fpath("contmodel/OrionSourceI_data_stretched_{0}.pdf".format(band)), bbox_inches='tight')
    fig5.savefig(paths.fpath("contmodel/OrionSourceI_data_stretched_{0}.svg".format(band)), bbox_inches='tight')




    fig6 = pl.figure(6)
    fig6.clf()
    ax1 = fig6.add_subplot(3,3,1)
    im = ax1.imshow(bestdiskmod_beam*jtok.value, interpolation='none',
                    origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm())
    vmin,vmax = im.norm.vmin, im.norm.vmax
    ax2 = fig6.add_subplot(3,3,2)
    im = ax2.imshow(bestdiskmod*jtok.value, interpolation='none',
                    origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    ax3 = fig6.add_subplot(3,3,3)
    im = ax3.imshow(bestdiskplussmearedsourcemod*jtok.value,
                    interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    cb = fig6.colorbar(mappable=im)
    cb.set_label("$T_B$ [K]")

    ax4 = fig6.add_subplot(3,3,4)
    im = ax4.imshow((data - bestdiskmod_beam)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    ax5 = fig6.add_subplot(3,3,5)
    im = ax5.imshow((data - bestdiskmod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    ax6 = fig6.add_subplot(3,3,6)
    im = ax6.imshow((data - bestdiskplussmearedsourcemod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    cb = fig6.colorbar(mappable=im)
    cb.set_label("$T_B$ [K]")

    ax7 = fig6.add_subplot(3,3,7)
    im = ax7.imshow((data - bestdiskmod_beam)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    vmin=-40, vmax=50,
                   )
    vmin,vmax = im.norm.vmin, im.norm.vmax
    ax8 = fig6.add_subplot(3,3,8)
    im = ax8.imshow((data - bestdiskmod_refit)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    vmin=vmin, vmax=vmax)
    ax9 = fig6.add_subplot(3,3,9)
    im = ax9.imshow((data - bestdiskplussmearedsourcemod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    vmin=vmin, vmax=vmax)
    cb = fig6.colorbar(mappable=im)
    cb.set_label("$T_B$ [K]")


    for ax in fig6.get_axes():
        ax.set_xticks([])
        ax.set_yticks([])

    pl.subplots_adjust(wspace=0, hspace=0.05)

    fig6.savefig(paths.fpath("contmodel/models_and_residuals_{0}.pdf".format(band)), bbox_inches='tight')


def crd_or_qty(x):
    try:
        return u.Quantity(x)
    except TypeError:
        return coordinates.SkyCoord(x)

freqs = list(fit_results.keys())
resultkeys = list(fit_results[freqs[0]].keys())
tabledata = [table.Column(data=u.Quantity(freqs), name='Frequency',)]
tabledata += [table.Column(data=crd_or_qty([fit_results[freq][key]
                                            for freq in freqs]),
                           name=key)
              for key in resultkeys]
fitted_coords = coordinates.SkyCoord([fit_results[freq]['Pt Position'] for freq in freqs])

tbl = table.Table(tabledata)
tbl['Pt Amp'] *= 1000
tbl['Pt Amp'].unit = u.mJy
tbl['ePt Amp'] *= 1000
tbl['ePt Amp'].unit = u.mJy
tbl['Total Flux'] *= 1000
tbl['Total Flux'].unit = u.mJy
tbl['Pt Width'].unit = u.au
tbl['Pt Flux'] *= 1000
tbl['Pt Flux'].unit = u.mJy

tbl['Disk PA'] = tbl['Disk PA'].to(u.deg)


formats = {'Pt Position': lambda x: "-" if np.isnan(x.ra) else "{0:0.4} {1:0.3}".format(x.ra.hms.s-14, x.dec.dms.s+30),
           'Pt RA': lambda x: "-" if np.isnan(x) else "{0:0.3}".format(x-14),
           'Pt Dec': lambda x: "-" if np.isnan(x) else "{0:0.3}".format(x+30),
           'ePt RA': lambda x: "-" if np.isnan(x) else "{0:0.4}".format(x),
           'ePt Dec': lambda x: "-" if np.isnan(x) else "{0:0.3}".format(x),
           'Pt Amp': lambda x: ('{0:0.2f}'.format(round_to_n(x,2))),
           'ePt Amp': lambda x: ('{0:0.2f}'.format(round_to_n(x,2))),
           'Pt Width': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk PA': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk FWHM': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk Radius': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Total Flux': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Pt Flux': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Pt \%': lambda x: strip_trailing_zeros('{0:0.5g}\%'.format(round_to_n(x,2)*100)),
          }


latexdict = latex_info.latexdict.copy()
latexdict['header_start'] = '\label{tab:continuum_fit_parameters}'
latexdict['caption'] = 'Continuum Fit Parameters'
latexdict['preamble'] = '\centering'
latexdict['tablefoot'] = ('\n\par The pointlike source '
                          'position is given as RA seconds and Dec arcseconds '
                          'offset from ICRS 5h35m14s -5d22m30s.  The error on '
                          'this position is 0.003s (RA) and 0.0003\\arcsec (Dec). '
                          'For the 7 mm data, the '
                          'position is left blank because we do not have '
                          'astrometric information for those data (they were '
                          'self-calibrated on a bright maser whose position '
                          'was not well-constrained).  '
                          'The disk FWHM'
                          ' is the vertical full-width half-maximum of the '
                          'fitted Gaussian profile.  '
                          'No formal parameter errors were measured for several of the '
                          '224.0 GHz fitted quantities because of a linear algebra '
                          'failure in the fitter; the errors are likely similar to '
                          'the 93.3 GHz fit errors. '
                          'The Pt Flux and Total Flux columns report the integrals '
                          'of the best-fit models.'
                         )

mask7mm = tbl['Frequency'] == 43.165
tbl['Pt Position'][mask7mm] = coordinates.SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg), frame='icrs')

tbl['Pt RA'] = table.Column(name='Pt RA', data=fitted_coords.ra.hms.s, unit=u.s)
tbl['Pt Dec'] = table.Column(name='Pt Dec', data=fitted_coords.dec.hms.s, unit=u.arcsec)
tbl['Pt RA'][mask7mm] = np.nan
tbl['Pt Dec'][mask7mm] = np.nan
tbl.remove_column('Pt Position')

tbl.sort('Frequency')
tbl.write('continuum_fit_parameters.txt', format='ascii.ecsv')

ntbl = table.Table.read('continuum_fit_parameters.txt', format='ascii.ecsv')
newformats = {k:v for k,v in formats.items()}


# these don't format well, so let's just skip them
ntbl.remove_column('ePt RA')
ntbl.remove_column('ePt Dec')
#ntbl.remove_column('ePt Amp')

for cn in ntbl.colnames:
    if 'e'+cn in ntbl.colnames:
        print("replacing {0}".format(cn))
        new_col = [("{0} $\pm$ {1}".format(formats[cn](val),
                                           (formats['e'+cn](err)
                                            if 'e'+cn in formats else
                                            formats[cn](err))
                                          )
                    if '-' != formats[cn](val) else '-') if err != 0.0
                   else formats[cn](val)
                   for val, err in zip(ntbl[cn], ntbl['e'+cn])]
        unit = ntbl[cn].unit
        ntbl.remove_column(cn)
        ntbl.remove_column('e'+cn)
        ntbl.add_column(table.Column(data=new_col, name=cn, unit=unit))
        newformats[cn] = lambda x: x


ntbl = ntbl['Frequency', 'Disk FWHM', 'Disk Radius', 'Disk PA', 'Pt RA', 'Pt Dec', 'Pt Amp', 'Pt Width', 'Pt Flux', 'Total Flux', 'Pt \%', ]
ntbl.sort('Frequency')
ntbl.write(paths.texpath('continuum_fit_parameters.tex'), format='ascii.latex',
           formats=newformats,
           latexdict=latexdict, overwrite=True)


with open(paths.texpath('continuum_beams.tex'), 'w') as fh:
    for band in beams:
        bm = beams[band]
        argdict = {'bandid': band.replace("3","three").replace("6", "six").replace("7", "seven"),
                   'major': strip_trailing_zeros('{0:0.5g}'.format(round_to_n(bm.major.to(u.arcsec).value,2))),
                   'minor': strip_trailing_zeros('{0:0.5g}'.format(round_to_n(bm.minor.to(u.arcsec).value,2))),
                   'pa': strip_trailing_zeros('{0:0.5g}'.format(round_to_n(bm.pa.to(u.deg).value,3))),
                  }
        fh.write("\\newcommand{{\\{bandid}maj}}{{{major}}}\n".format(**argdict))
        fh.write("\\newcommand{{\\{bandid}min}}{{{minor}}}\n".format(**argdict))
        fh.write("\\newcommand{{\\{bandid}pa}}{{{pa}}}\n".format(**argdict))

print()
print("Measured point source offsets: ")
print(ptsrc_diskcen_offsets)
for k in ['B6','B3']:
    print("{0} fitted offset and positional error".format(k))
    print((ptsrc_diskcen_offsets[k][0]*d_orion).to(u.au, u.dimensionless_angles()),
          (ptsrc_diskcen_offsets[k][2]*d_orion).to(u.au, u.dimensionless_angles()))

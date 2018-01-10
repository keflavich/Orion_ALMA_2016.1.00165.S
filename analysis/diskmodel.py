# attempt to model Source I disk
import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from astropy import convolution
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

fit_results = {}

beams = {}
 
for fn, freq, band in [#('Orion_SourceI_B6_continuum_r-2_longbaselines_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                       #('Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                       ('Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4_SourceIcutout.image.tt0.pbcor.fits', 224.0*u.GHz, 'B6'),
                       ('Orion_SourceI_B3_continuum_r-2.clean0.25mJy_SourceIcutout.image.tt0.pbcor.fits', 93.3*u.GHz, 'B3'),
                      ]:

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
    jtok = observed_beam.jtok(freq)

    def model(x1, x2, y1, y2, scale, kernelmajor=None, kernelminor=None, kernelpa=None,
              ptsrcx=None, ptsrcy=None, ptsrcamp=None):
        #if any(ii < 0 for ii in (x1,x2,y1,y2)):
        #    return 1e5
        #if (x1 > data.shape[1]-1) or (x2 > data.shape[1]-1) or (y1 > data.shape[0]-1) or (y2 > data.shape[0]-1):
        #    return 1e5
        if x2 > data.shape[1] - 1:
            x2 = data.shape[1] -1
        if y2 > data.shape[0] - 1:
            y2 = data.shape[0] -1

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

        beam_kernel = beam.as_kernel(pixscale)
        beam_amp = beam_kernel.array.max()

        diskmod = convolution.convolve_fft(disk, beam_kernel) / beam_amp

        if ptsrcamp is not None:
            ptsrcmod = shift.shift2d(observed_beam.as_kernel(pixscale,
                                                             x_size=data.shape[1],
                                                             y_size=data.shape[0]),
                                     ptsrcx-data.shape[0]/2, ptsrcy-data.shape[1]/2) / beam_amp * ptsrcamp
            diskmod += ptsrcmod


        return diskmod

    def residual(pars):
        mod = model(**pars)
        return (data - mod)

    diskmod = model(x1,x2,y1,y2,data.max())

    parameters = lmfit.Parameters()
    parameters.add('x1', value=x1)
    parameters.add('x2', value=x2)
    parameters.add('y1', value=y1)
    parameters.add('y2', value=y2)
    parameters.add('scale', value=data.max())
    result = lmfit.minimize(residual, parameters, epsfcn=0.05)
    print("Basic fit parameters (linear model):")
    print(result.params)
    print()

    bestdiskmod_beam = model(**result.params)

    # Create a "beam" that is really the vertical x horizontal scale height
    # to be convolved with the observed beam
    # old versions for QA2 cutout parameters.add('kernelmajor', value=0.064)
    # old versions for QA2 cutout parameters.add('kernelminor', value=0.043)
    parameters.add('kernelmajor', value=0.054)
    parameters.add('kernelminor', value=0.033)
    # Fix the position angle such that one direction of the resulting kernel will
    # directly be a Gaussian scale height
    measured_positionangle = 37
    parameters.add('kernelpa', value=measured_positionangle+90, vary=False)
    result2 = lmfit.minimize(residual, parameters, epsfcn=0.1)
    print("Smoothed linear fit parameters:")
    print(result2.params)
    print()

    bestdiskmod = model(**result2.params)

    # from the first round of fitting, there is a residual source at this position
    ptsrc = coordinates.SkyCoord(83.81049240934931*u.deg, -5.375170355557261*u.deg, frame='icrs')

    ptsrcx, ptsrcy = mywcs.wcs_world2pix(ptsrc.ra, ptsrc.dec, 0)
    ptsrc_amp_value = data[int(ptsrcy), int(ptsrcx)]
    print("Point source amplitude data value: {0}".format(ptsrc_amp_value))
    print()

    parameters.add('ptsrcx', value=ptsrcx, min=ptsrcx-5, max=ptsrcx+5)
    parameters.add('ptsrcy', value=ptsrcy, min=ptsrcy-5, max=ptsrcy+5)
    parameters.add('ptsrcamp', value=0.004, min=0.001, max=0.6)
    result3 = lmfit.minimize(residual, parameters, epsfcn=0.1)
    print("Smoothed linear fit parameters with point source:")
    print(result3.params)
    print()

    bestdiskplussourcemod = model(**result3.params)

    ptsrc_ra, ptsrc_dec = mywcs.wcs_pix2world(result3.params['ptsrcx'], result3.params['ptsrcy'], 0)
    fitted_ptsrc = coordinates.SkyCoord(ptsrc_ra*u.deg, ptsrc_dec*u.deg, frame=mywcs.wcs.radesys.lower())
    print("Fitted point source location = {0} {1}".format(fitted_ptsrc.to_string('hmsdms'), fitted_ptsrc.frame.name))
    print("Fitted point source amplitude: {0}".format(result3.params['ptsrcamp']))

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

    posang = np.arctan2(result3.params['y2']-result3.params['y1'],
                        result3.params['x2']-result3.params['x1'])*u.rad - 90*u.deg
    print("posang={0}".format(posang.to(u.deg)))

    fitted_beam = radio_beam.Beam(result2.params['kernelmajor']*u.arcsec,
                                  result2.params['kernelminor']*u.arcsec,
                                  result2.params['kernelpa']*u.deg,)
    # NOTE: the fitted beam *is* the source size after a revision to the model
    # in which the input beam is convolved with the observed beam
    #source_size = fitted_beam.deconvolve(observed_beam)
    source_size = fitted_beam
    print("Fitted source (disk vertical scale) size: {0}".format(fitted_beam.__repr__()))
    print("Real source (disk vertical scale) size: {0}".format(source_size.__repr__()))
    scaleheight = (fitted_beam.major*d_orion).to(u.au, u.dimensionless_angles())
    print("Scale height: {0}".format(scaleheight))

    length_as = (((result2.params['x2'] - result2.params['x1'])**2 +
                  (result2.params['y2']-result2.params['y1'])**2)**0.5 * pixscale).to(u.arcsec)
    length_au = (length_as * d_orion).to(u.au, u.dimensionless_angles())

    print("Length in arcsec: {0:0.3g}  in AU: {1:0.3g}  or radius {2:0.3g}"
          .format(length_as, length_au, length_au/2))
    print()
    print()

    ppbeam = (observed_beam.sr/pixscale**2).decompose()

    fit_results[freq] = {
                         'Disk Scaleheight': scaleheight,
                         'Disk Radius': length_au/2,
                         'Disk PA': posang,
                         'Ptsrc Position': fitted_ptsrc,
                         'Ptsrc Amp': result3.params['ptsrcamp'].value,
                         'Total Flux': bestdiskplussourcemod.sum() / ppbeam,
                         'Ptsrc Frac': result3.params['ptsrcamp'].value / (bestdiskplussourcemod.sum() / ppbeam),
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
    pl.imshow(data_K, interpolation='none', origin='lower', cmap='viridis')
    pl.colorbar()



    # publication figures
    fig5 = pl.figure(5)
    fig5.clf()
    ax = fig5.gca()
    im0 = ax.imshow(data*jtok.value, interpolation='none', origin='lower', cmap='viridis')
    im = ax.imshow(data*1e3, interpolation='none', origin='lower', cmap='viridis')
    cb2 = fig5.colorbar(mappable=im0)
    cb2.set_label('$T_B$ [K]')
    cb = fig5.colorbar(mappable=im)
    cb.set_label('$S_{1.3 mm}$ [mJy beam$^{-1}$]')
    ax.set_xticks([])
    ax.set_yticks([])

    beam_ellipse = observed_beam.ellipse_to_plot(0.9*data.shape[1], 0.1*data.shape[0], pixscale)
    beam_ellipse.set_edgecolor('w')
    beam_ellipse.set_facecolor('w')
    ax.add_patch(beam_ellipse)
    fig5.savefig(paths.fpath("contmodel/OrionSourceI_data_{0}.pdf".format(band)), bbox_inches='tight')

    fig5.clf()
    ax = fig5.gca()
    im0 = ax.imshow(data*jtok.value, interpolation='none', origin='lower', cmap='viridis', vmin=-5, vmax=40)
    im = ax.imshow(data*1e3, interpolation='none', origin='lower', cmap='viridis', vmin=-5/jtok.value*1e3, vmax=40/jtok.value*1e3)
    cb2 = fig5.colorbar(mappable=im0)
    cb2.set_label('$T_B$ [K]')
    cb = fig5.colorbar(mappable=im)
    cb.set_label('$S_{1.3 mm}$ [mJy beam$^{-1}$]')
    ax.set_xticks([])
    ax.set_yticks([])

    # need duplicate code here (grr) because matplotlib refuses to reuse Aritst objects
    beam_ellipse = observed_beam.ellipse_to_plot(0.9*data.shape[1], 0.1*data.shape[0], pixscale)
    beam_ellipse.set_edgecolor('w')
    beam_ellipse.set_facecolor('w')
    ax.add_patch(beam_ellipse)
    fig5.savefig(paths.fpath("contmodel/OrionSourceI_data_stretched_{0}.pdf".format(band)), bbox_inches='tight')




    fig6 = pl.figure(6)
    fig6.clf()
    ax1 = fig6.add_subplot(3,3,1)
    im = ax1.imshow(bestdiskmod_beam*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm())
    vmin,vmax = im.norm.vmin, im.norm.vmax
    ax2 = fig6.add_subplot(3,3,2)
    im = ax2.imshow(bestdiskmod*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    norm=asinh_norm.AsinhNorm(),
                    vmin=vmin, vmax=vmax)
    ax3 = fig6.add_subplot(3,3,3)
    im = ax3.imshow(bestdiskplussourcemod*jtok.value, interpolation='none', origin='lower', cmap='viridis',
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
    im = ax6.imshow((data - bestdiskplussourcemod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
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
    im = ax8.imshow((data - bestdiskmod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
                    vmin=vmin, vmax=vmax)
    ax9 = fig6.add_subplot(3,3,9)
    im = ax9.imshow((data - bestdiskplussourcemod)*jtok.value, interpolation='none', origin='lower', cmap='viridis',
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
tbl = table.Table(tabledata)
tbl['Ptsrc Amp'] *= 1000
tbl['Ptsrc Amp'].unit = u.mJy
tbl['Total Flux'] *= 1000
tbl['Total Flux'].unit = u.mJy

tbl['Disk PA'] = tbl['Disk PA'].to(u.deg)


formats = {'Ptsrc Position': lambda x: x.to_string('hmsdms'),
           'Ptsrc Amp': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk PA': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk Scaleheight': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Disk Radius': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Total Flux': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           'Ptsrc Frac': lambda x: strip_trailing_zeros('{0:0.5g}'.format(round_to_n(x,2))),
          }


latexdict = latex_info.latexdict.copy()
latexdict['header_start'] = '\label{tab:continuum_fit_parameters}'
latexdict['caption'] = 'Continuum Fit Parameters'
latexdict['preamble'] = '\centering'
latexdict['tablefoot'] = ('\n\par The centroid coordinate of the point source '
                          'is specified in ICRS coordinates.')



tbl.write(paths.texpath('continuum_fit_parameters.tex'), format='ascii.latex',
          formats=formats,
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

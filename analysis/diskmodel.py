# attempt to model Source I disk
import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from astropy import convolution
import radio_beam
import lmfit
from constants import d_orion
from image_registration.fft_tools import shift

fh = fits.open('/Users/adam/work/orion/alma_lb/FITS/uid_A001_X88e_X1d3_calibrated_final_cont.pbcor.fits')

cropslice_x = slice(3510,3610)
cropslice_y = slice(3520,3620)
data = fh[0].data.squeeze()[cropslice_y, cropslice_x]
mywcs = wcs.WCS(fh[0].header).celestial[cropslice_y, cropslice_x]
pixscale = wcs.utils.proj_plane_pixel_area(mywcs)**0.5 * u.deg

diskends = coordinates.SkyCoord(['5:35:14.5232 -5:22:30.73',
                                 '5:35:14.5132 -5:22:30.54'],
                                frame='fk5', unit=(u.hour, u.deg))

diskends_pix = np.array(mywcs.wcs_world2pix(diskends.ra.deg, diskends.dec.deg, 0))
(x1,x2),(y1,y2) = diskends_pix

observed_beam = radio_beam.Beam.from_fits_header(fh[0].header)

data_K = (data*u.Jy).to(u.K, observed_beam.jtok_equiv(fh[0].header['CRVAL3']*u.Hz))

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
result = lmfit.minimize(residual, parameters, epsfcn=0.1)
print(result.params)

bestdiskmod_beam = model(**result.params)

# Create a "beam" that is really the vertical x horizontal scale height
# to be convolved with the observed beam
parameters.add('kernelmajor', value=0.064)
parameters.add('kernelminor', value=0.043)
# Fix the position angle such that one direction of the resulting kernel will
# directly be a Gaussian scale height
# (37 deg is the measured PA)
parameters.add('kernelpa', value=37+90, vary=False)
result2 = lmfit.minimize(residual, parameters, epsfcn=0.1)
print(result2.params)

bestdiskmod = model(**result2.params)

# from the first round of fitting, there is a residual source at this position
ptsrc = coordinates.SkyCoord(83.81049240934931*u.deg, -5.375170355557261*u.deg, frame='icrs')

ptsrcx, ptsrcy = mywcs.wcs_world2pix(ptsrc.ra, ptsrc.dec, 0)

parameters.add('ptsrcx', value=ptsrcx, min=ptsrcx-5, max=ptsrcx+5)
parameters.add('ptsrcy', value=ptsrcy, min=ptsrcy-5, max=ptsrcy+5)
parameters.add('ptsrcamp', value=0.004, min=0.002, max=0.6)
result3 = lmfit.minimize(residual, parameters, epsfcn=0.1)
print(result3.params)

bestdiskplussourcemod = model(**result3.params)

ptsrc_ra, ptsrc_dec = mywcs.wcs_pix2world(result3.params['ptsrcx'], result3.params['ptsrcy'], 0)
fitted_ptsrc = coordinates.SkyCoord(ptsrc_ra*u.deg, ptsrc_dec*u.deg, frame=mywcs.wcs.radesys.lower())
print("Fitted point source location = {0}".format(fitted_ptsrc.to_string('hmsdms')))

import pylab as pl
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

pl.savefig("SourceI_Disk_model.png")

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

pl.savefig("SourceI_Disk_model_bigbeam.png")

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

pl.savefig("SourceI_Disk_model_bigbeam_withptsrc.png")

pl.figure(4)
pl.clf()
pl.imshow(data_K, interpolation='none', origin='lower', cmap='viridis')
pl.colorbar()


posang = np.arctan2(result.params['x2']-result.params['x1'],
                    result.params['y2']-result.params['y1'])*u.rad
print("posang={0}".format(posang.to(u.deg)))

fitted_beam = radio_beam.Beam(result2.params['kernelmajor']*u.arcsec,
                              result2.params['kernelminor']*u.arcsec,
                              result2.params['kernelpa']*u.deg,)
# NOTE: the fitted beam *is* the source size after a revision to the model
# in which the input beam is convolved with the observed beam
#source_size = fitted_beam.deconvolve(observed_beam)
source_size = fitted_beam
print("Fitted source size: {0}".format(fitted_beam.__repr__()))
print("Real source size: {0}".format(source_size.__repr__()))
scaleheight = (fitted_beam.major*d_orion).to(u.au, u.dimensionless_angles())
print("Scale height: {0}".format(scaleheight))

length_as = (((result2.params['x2'] - result2.params['x1'])**2 +
              (result2.params['y2']-result2.params['y1'])**2)**0.5 * pixscale).to(u.arcsec)
length_au = (length_as * d_orion).to(u.au, u.dimensionless_angles())

print("Length in arcsec: {0:0.3g}  in AU: {1:0.3g}  or radius {2:0.3g}"
      .format(length_as, length_au, length_au/2))

# attempt to model Source I disk
import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from astropy import convolution
import radio_beam
import lmfit

fh = fits.open('/Users/adam/work/orion/alma/FITS/longbaseline/uid_A001_X88e_X1d3_calibrated_final_cont.pbcor.fits')

cropslice = slice(3510,3630)
data = fh[0].data.squeeze()[cropslice, cropslice]
mywcs = wcs.WCS(fh[0].header).celestial[cropslice,cropslice]
pixscale = wcs.utils.proj_plane_pixel_area(mywcs)**0.5 * u.deg

diskends = coordinates.SkyCoord(['5:35:14.5232 -5:22:30.73',
                                 '5:35:14.5132 -5:22:30.54'],
                                frame='fk5', unit=(u.hour, u.deg))

diskends_pix = np.array(mywcs.wcs_world2pix(diskends.ra.deg, diskends.dec.deg, 0))
(x1,x2),(y1,y2) = diskends_pix

observed_beam = radio_beam.Beam.from_fits_header(fh[0].header)

def model(x1, x2, y1, y2, scale, kernelmajor=None, kernelminor=None, kernelpa=None):
    def line_y(x):
        # y = m x + b
        m = (y2-y1)/(x2-x1)
        b = y1 - x1*m
        return m*x+b

    xx = np.linspace(x1, x2, 1000)
    yy = line_y(xx)

    disk = np.zeros_like(data, dtype='float')
    disk[(np.round(yy).astype('int')), (np.round(xx).astype('int'))] = 1.

    if kernelmajor is None:
        beam = observed_beam
    else:
        beam = radio_beam.Beam(kernelmajor*u.arcsec, kernelminor*u.arcsec, kernelpa*u.deg)
    diskmod = convolution.convolve_fft(disk, beam.as_kernel(pixscale))

    return diskmod/diskmod.max() * scale

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

bestdiskmod_beam = model(**result.params)

parameters.add('kernelmajor', value=0.052)
parameters.add('kernelminor', value=0.03)
parameters.add('kernelpa', value=-77)
result2 = lmfit.minimize(residual, parameters, epsfcn=0.1)

bestdiskmod = model(**result2.params)

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



posang = np.arctan2(result.params['x2']-result.params['x1'],
                    result.params['y2']-result.params['y1'])*u.rad
print("posang={0}".format(posang.to(u.deg)))

fitted_beam = radio_beam.Beam(result2.params['kernelmajor']*u.arcsec,
                              result2.params['kernelminor']*u.arcsec,
                              result2.params['kernelpa']*u.deg,)
source_size = fitted_beam.deconvolve(observed_beam)
print("Real source size: {0}".format(source_size.__repr__()))

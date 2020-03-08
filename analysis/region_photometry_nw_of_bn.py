import numpy as np
import files
import regions
import paths
import constants
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import radio_beam


ring_ap_fn = paths.rpath('nw_of_bn_ring.reg')
regfs = {x: ring_ap_fn for x in (3,6,7)}
regfs[9] = paths.rpath('nw_of_bn_ring_660GHzAp.reg')
regfs[8] = paths.rpath('nw_of_bn_ring_660GHzAp.reg')
regs = {x: regions.read_ds9(regfs[x]) for x in regfs}

photdata = {}
photpeak = {}
photstd = {}
photsum = {}
beams = {}
for fn,band in zip([files.b3_hires_cont, files.b6_hires_cont, files.b7_hires_cont, files.b8_cont, files.b9_cont],
                   (3,6,7,8,9),
                  ):
    fh = fits.open(paths.dpath(fn))
    ww = wcs.WCS(fh[0].header).celestial
    preg = regs[band][0].to_pixel(ww)
    msk = preg.to_mask()
    data = fh[0].data.squeeze()
    mdata = msk.multiply(data)
    mdata = mdata[mdata != 0]
    photsum[band] = mdata.sum()
    photpeak[band] = mdata.max()
    photstd[band] = mdata.std()
    beams[band] = radio_beam.Beam.from_fits_header(fh[0].header)
    photdata[band] = photsum[band] * u.Jy / (beams[band].sr / u.Quantity(wcs.utils.proj_plane_pixel_area(ww), u.deg**2)).decompose()

print(photdata)

import pylab as pl
pl.clf()
frqs = u.Quantity([constants.central_freqs[f'B{bnd}'] for bnd in (3,6,7,9)])
pl.loglog(frqs,
          u.Quantity([photdata[bnd] for bnd in (3,6,7,9)]),
          's')
pl.loglog(constants.central_freqs['B8'],
          u.Quantity([photpeak[bnd] for bnd in (8,)]),
          'v')
pl.loglog(constants.central_freqs['B8'],
          u.Quantity([photstd[bnd] for bnd in (8,)]),
          'v')

linfrqs = np.linspace(frqs[0], frqs[-1])
pl.plot(linfrqs, (linfrqs / frqs[0])**2 * photdata[3], 'k:')
#pl.plot(linfrqs, (linfrqs / frqs[1])**2 * photdata[6], 'k:')
#pl.plot(linfrqs, (linfrqs / frqs[2])**2 * photdata[7], 'k:')
#pl.plot(linfrqs, (linfrqs / frqs[0])**3.5 * photdata[3], 'r--')
#pl.plot(linfrqs, (linfrqs / frqs[1])**3.5 * photdata[6], 'r--')
pl.plot(linfrqs, (linfrqs / frqs[2])**3.5 * photdata[7], 'r--')
pl.xlabel("Frequency (GHz)")
pl.ylabel("Flux Density (Jy)")

pl.savefig(paths.fpath("nw_of_bn_sed.png"))

import dust_emissivity
print("B3 mass: ",dust_emissivity.dust.massofsnu(constants.central_freqs['B3'], photdata[3], distance=400*u.pc, beta=1.5))
print("B6 mass: ",dust_emissivity.dust.massofsnu(constants.central_freqs['B6'], photdata[6], distance=400*u.pc, beta=1.5))
print("B7 mass: ",dust_emissivity.dust.massofsnu(constants.central_freqs['B7'], photdata[7], distance=400*u.pc, beta=1.5))
print("B9 mass: ",dust_emissivity.dust.massofsnu(constants.central_freqs['B9'], photdata[9], distance=400*u.pc, beta=1.5))

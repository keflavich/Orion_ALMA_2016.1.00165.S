import os
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import log
import shutil
from astropy.io import fits
import pyregion
import paths
import pylab as pl
import pyspeckit



fn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout.fits'
outfn = '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout_K.fits'
outfn_medsub= '/Volumes/external/orion/full_OrionSourceI_B6_spw0_lines_cutout_K_medsub.fits'

def make_kelvin_cube():
    cube = SpectralCube.read(fn)

    # this copy step is necessary to allocate memory for the output
    shutil.copy(fn, outfn)
    outfh = fits.open(outfn, mode='update')

    jtok_factors = cube.jtok_factors()
    for index,(slice,factor) in enumerate(zip(cube,jtok_factors)):
        outfh[0].data[index] = slice * factor
        outfh.flush()

    outfh[0].header['BUNIT'] = 'K'
    outfh.flush()
    outfh.close()

    #cube = (SpectralCube.read(fn)[:,515:721,550:714].mask_out_bad_beams(5))



def make_medsub_cube():
    log.info("Calculating median")
    cube = SpectralCube.read(outfn)
    cube = cube.mask_out_bad_beams(0.01)
    med = cube.percentile(25,axis=0)

    shutil.copy(outfn, outfn_medsub)

    outfh = fits.open(outfn_medsub, mode='update')

    for index,slice in enumerate(cube):
        outfh[0].data[index] = slice - med
        outfh.flush()

    outfh.flush()
    outfh.close()

cube = SpectralCube.read(outfn)
cofn = '/Volumes/external/orion/12CO2-1_OrionSourceI.fits'
if not os.path.exists(cofn):
    log.info("Writing CO cube")
    cocube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=230.538*u.GHz).spectral_slab(-100*u.km/u.s, 120*u.km/u.s)
    cocube.write(cofn, overwrite=True)

medsub = SpectralCube.read(outfn_medsub)
log.info("Computing spectral mean")
medsub_mean = medsub.mean(axis=(1,2), how='slice', progressbar=True)
medsub_mean.write('/Volumes/external/orion/OrionSourceI_medsub_mean_spectrum.fits', overwrite=True)
medsub_max = medsub.max(axis=(1,2), how='slice', progressbar=True)
medsub_max.write('/Volumes/external/orion/OrionSourceI_medsub_max_spectrum.fits', overwrite=True)


reg = pyregion.open(paths.rpath('sourceI_innerdisk_box.reg'))
subcube = cube.subcube_from_ds9region(reg)
sourceIdiskspec = subcube.mean(axis=(1,2))
sourceIdiskspec.write('/Volumes/external/orion/OrionSourceI_innerdisk_mean_spectrum.fits', overwrite=True)

sp1 = pyspeckit.Spectrum('/Volumes/external/orion/OrionSourceI_medsub_mean_spectrum.fits')
sp1.plotter()
sp1.plotter.savefig(paths.fpath('OrionSourceI_medsub_mean_spectrum.png'))

sp2 = pyspeckit.Spectrum('/Volumes/external/orion/OrionSourceI_medsub_max_spectrum.fits')
sp2.plotter()
sp2.plotter.savefig(paths.fpath('OrionSourceI_medsub_max_spectrum.png'))

sp3 = pyspeckit.Spectrum('/Volumes/external/orion/OrionSourceI_innerdisk_mean_spectrum.fits')
sp3.plotter()
sp3.plotter.savefig(paths.fpath('OrionSourceI_innerdisk_mean_spectrum.png'))

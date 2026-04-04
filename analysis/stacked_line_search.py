"""
Script to create a velocity map (moment 1 map) from a known "good" line, in
this case one of the brighter NaCl transitions, then use that velocity map to
shift-and-stack all of the spectra in the cube.

The main functionality is spectral-cube's stacking function:
https://github.com/radio-astro-tools/spectral-cube/blob/master/spectral_cube/analysis_utilities.py#L136

"""
import numpy as np
import os
import spectral_cube.analysis_utilities
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.io import fits
import pylab as pl
import regions
import reproject

# the 'paths' module specifies paths to the filenames.  It is specific to the
# Orion ALMA project; you need to either modify it or remove these imports and
# replace calls to the 'paths' functions below
import paths
from paths import fcp


basedir = '/orange/adamginsburg/orion'
project_2025_1_00236 = [
    ('B4', 25, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da5/group.uid___A001_X3833_X4da6/member.uid___A001_X3833_X4da7/product/member.uid___A001_X3833_X4da7.Orion_SrcI_sci.spw25.cube.I.selfcal.pbcor.fits'),
    ('B4', 27, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da5/group.uid___A001_X3833_X4da6/member.uid___A001_X3833_X4da7/product/member.uid___A001_X3833_X4da7.Orion_SrcI_sci.spw27.cube.I.selfcal.pbcor.fits'),
    ('B4', 29, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da5/group.uid___A001_X3833_X4da6/member.uid___A001_X3833_X4da7/product/member.uid___A001_X3833_X4da7.Orion_SrcI_sci.spw29.cube.I.selfcal.pbcor.fits'),
    ('B4', 31, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da5/group.uid___A001_X3833_X4da6/member.uid___A001_X3833_X4da7/product/member.uid___A001_X3833_X4da7.Orion_SrcI_sci.spw31.cube.I.selfcal.pbcor.fits'),
    ('B4', 33, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da5/group.uid___A001_X3833_X4da6/member.uid___A001_X3833_X4da7/product/member.uid___A001_X3833_X4da7.Orion_SrcI_sci.spw33.cube.I.selfcal.pbcor.fits'),
    ('B6high', 25, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da9/group.uid___A001_X3833_X4daa/member.uid___A001_X3833_X4dab/product/member.uid___A001_X3833_X4dab.Orion_SrcI_sci.spw25.cube.I.pbcor.fits'),
    ('B6high', 27, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da9/group.uid___A001_X3833_X4daa/member.uid___A001_X3833_X4dab/product/member.uid___A001_X3833_X4dab.Orion_SrcI_sci.spw27.cube.I.pbcor.fits'),
    ('B6high', 29, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da9/group.uid___A001_X3833_X4daa/member.uid___A001_X3833_X4dab/product/member.uid___A001_X3833_X4dab.Orion_SrcI_sci.spw29.cube.I.pbcor.fits'),
    ('B6high', 31, f'{basedir}/2025.1.00236.S/2025.1.00236.S/science_goal.uid___A001_X3833_X4da9/group.uid___A001_X3833_X4daa/member.uid___A001_X3833_X4dab/product/member.uid___A001_X3833_X4dab.Orion_SrcI_sci.spw31.cube.I.pbcor.fits'),
]


# step 1: create a velocity map

vmap_name = paths.dpath('disk_velocity_map.fits')
if not os.path.exists(vmap_name):
    cube = SpectralCube.read(paths.dpath('cubes/OrionSourceI_Unknown_4_robust0.5.maskedclarkclean10000_medsub_K.fits'))
    cube = SpectralCube.read(paths.dpath('cubes/OrionSourceI_Unknown_4_robust0.5maskedclarkclean10000_medsub_K.fits'))
    m1 = cube.moment1()
    m0 = cube.moment0()
    mask = m0.value > 300

    vmap = m1
    vmap[~mask] = np.nan

    r =regions.Regions.read(paths.rpath('sourceI_enclosing_ellipse.reg'))[0]
    rp = r.to_pixel(vmap.wcs)
    mask = rp.to_mask()

    vmap_ = np.empty(vmap.shape)*np.nan
    vmap_[mask.bbox.slices] = vmap[mask.bbox.slices].value * mask.data
    hdu = vmap.hdu
    hdu.data = vmap_
    hdu.writeto(vmap_name, overwrite=True)
else:
    hdu = fits.open(vmap_name)[0]
vmap = spectral_cube.lower_dimensional_structures.Projection.from_hdu(hdu)


# step 2: stack

cubes_to_process = []

for band in ('B3', 'B6', 'B7'):
    for spw in (0,1,2,3):
        for robust in (-2, 0.5, 2):

            suffix = '.lb' if band == 'B7' else ''

            fn = None
            templates = [
                'OrionSourceI_only.{1}{3}.robust{2}.spw{0}.maskedclarkclean10000_medsub.image.pbcor.cb.K.fits',
                'OrionSourceI_only.{1}{3}.robust{2}.spw{0}.clarkclean10000_medsub.image.pbcor.cb.K.fits',
                'OrionSourceI_only.{1}{3}.robust{2}.spw{0}.maskedclarkclean10000.image.pbcor.cb.K.fits',
                'OrionSourceI_only.{1}{3}.robust{2}.spw{0}.clarkclean10000.image.pbcor.cb.K.fits',
            ]
            for template in templates:
                candidate = fcp(template.format(spw, band, robust, suffix))
                external_candidate = candidate.replace('/imaging/', '/external/')
                if os.path.exists(candidate):
                    fn = candidate
                    break
                if os.path.exists(external_candidate):
                    fn = external_candidate
                    break
            if fn is None:
                print(f"Skipping missing cube for band={band} spw={spw} robust={robust}")
                continue

            outname = 'OrionSourceI_{1}{3}_spw{0}_robust{2}'.format(spw, band, robust, suffix)
            cubes_to_process.append((fn, outname))

for band, spw, fn in project_2025_1_00236:
    outname = 'OrionSourceI_2025.1.00236.S_{0}_spw{1}'.format(band, spw)
    cubes_to_process.append((fn, outname))

for fn, outname in cubes_to_process:
    fullcube = SpectralCube.read(fn, use_dask=True)
    print(fn, fullcube.spectral_extrema)

    # convert the cube to velocity units with an arbitrary reference point
    # (this step assumes the cube is in frequency or wavelength; if the
    # cube is not, it should be skipped)
    fullcube = fullcube.with_spectral_unit(u.km/u.s,
                                           velocity_convention='radio',
                                           rest_value=fullcube.spectral_axis.mean())

    # mask out super bright SiO masers; threshold depends on cube units
    if fullcube.unit.is_equivalent(u.Jy/u.beam):
        fullcube = fullcube.with_mask(fullcube < 0.5*u.Jy/u.beam)
    elif fullcube.unit.is_equivalent(u.K):
        fullcube = fullcube.with_mask(fullcube < 500*u.K)

    # reproject the velocity map into the cube's coordinate system
    vmap_proj,_ = reproject.reproject_interp(vmap.hdu,
                                             fullcube.wcs.celestial,
                                             shape_out=fullcube.shape[1:])
    vmap_proj = u.Quantity(vmap_proj, u.km/u.s)

    # perform the stacking!
    stack = spectral_cube.analysis_utilities.stack_spectra(fullcube, vmap_proj,
                                                           v0=0.0*u.km/u.s)
    fstack = stack.with_spectral_unit(u.GHz)

    fstack.write(paths.dpath(f'stacked_spectra/{outname}_K.fits'),
                 overwrite=True)

    pl.clf()
    fstack.quicklook(filename=paths.fpath(f'stacked_spectra/{outname}.pdf'))

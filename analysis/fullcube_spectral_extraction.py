import os
from astropy import log
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube
import pyregion
import radio_beam
import socket

if 'nmpost' in socket.gethostname():
    dpath = lambda x: os.path.join("/lustre/aginsbur/orion/2016.1.00165.S/imaging",x)
    rpath = lambda x: os.path.join("/lustre/aginsbur/orion/2016.1.00165.S/regions",x)
    spath = lambda x: os.path.join("/lustre/aginsbur/orion/2016.1.00165.S/imaging/spectra",x)
elif 'rng9000' in socket.gethostname():
    dpath = lambda x: os.path.join("/media/Seagate Expansion Drive-1/orion/2016.1.00165.S/imaging",x)
    rpath = lambda x: os.path.join("/lustre/aginsbur/orion/2016.1.00165.S/regions",x)
    spath = lambda x: os.path.join("/lustre/aginsbur/orion/2016.1.00165.S/imaging/spectra",x)
else:
    raise ValueError("No match to socket hostname {0}.".format(socket.gethostname()))

mergecubes = [
    'full_OrionSourceI_B3_spw0_lines.fits',
    'full_OrionSourceI_B3_spw3_lines.fits',
    'full_OrionSourceI_B3_spw2_lines.fits',
    'full_OrionSourceI_B3_spw1_lines.fits',
    'full_OrionSourceI_B7_maskedclean_spw1_lines.fits',
    'full_OrionSourceI_B7_maskedclean_spw0_lines.fits',
    'full_OrionSourceI_B7_maskedclean_spw2_lines.fits',
    'full_OrionSourceI_B7_maskedclean_spw3_lines.fits',
    'full_OrionSourceI_B6_maskedclean_spw3_lines.fits',
    'full_OrionSourceI_B6_maskedclean_spw2_lines.fits',
    'full_OrionSourceI_B6_maskedclean_spw1_lines.fits',
    'full_OrionSourceI_B6_maskedclean_spw0_lines.fits',
]


regions = (
           pyregion.open(rpath('apertures_B6_ref.reg'))
           #pyregion.open(rpath('tc_continuum_core_extraction_regions.reg')) +
           #pyregion.open(rpath('ionizationfront_circle.reg')) +
           #pyregion.open(rpath('extraction_regions_n_and_m.reg')) +
           #pyregion.open(rpath('ch3cn_large_cores.reg'))
          )

for cubename in mergecubes:
    print("Loading cube {0}".format(dpath(cubename)))
    cube = SpectralCube.read(dpath(cubename)).mask_out_bad_beams(threshold=0.1)

    for reg in regions:

        name = reg.attr[1]['text']
        fname = name.replace(" ","_").lower()

        suffix = os.path.splitext(cubename)[0]
        if os.path.exists(spath("{1}_{0}.fits".format(suffix,fname))):
            print("Skipping {0} {1} because it exists".format(suffix, fname))
            continue

        CL = reg.coord_list
        if reg.name == 'circle':
            radius = CL[2]*3600
            reg = pyregion.ShapeList([reg])
        elif reg.name == 'ellipse':
            radius = CL[2]*3600 # just use major
            reg = pyregion.ShapeList([reg])
        else:
            radius = 0.5
            reg = pyregion.parse("fk5; circle({0},{1},0.5\")"
                                 .format(CL[0], CL[1]))

        try:
            scube = cube.subcube_from_ds9region(reg)
        except ValueError as ex:
            print("Skipping {0} because {1}".format(name, ex))
            continue
        print(cube)
        log.info("Source name: {0}  filename: {1}".format(name,fname))
        print(scube)
        spsum = scube.sum(axis=(1,2))
        assert np.any(np.isfinite(spsum))
        spnpix = np.count_nonzero(np.isfinite(scube[1000,:,:]))
        assert spnpix > 0
        spectrum = spsum / spnpix
        # I think this is a hack left over from old versions of SpectralCube
        spsum.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in scube.beams]),
                                             minor=np.nanmedian([bm.minor.to(u.deg).value for bm in scube.beams]),
                                             pa=np.nanmedian([bm.pa.to(u.deg).value for bm in scube.beams]),)

        hdu = spsum.hdu
        hdu.data = spectrum.value
        pixel_scale = np.abs(cube.wcs.celestial.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        hdu.header['PPBEAM'] = (spsum.meta['beam'].sr / pixel_scale**2).decompose().value

        hdu.header['OBJECT'] = name

        hdu.writeto(spath("{1}_{0}.fits".format(suffix,fname)), clobber=True)
        print(spath("{1}_{0}.fits".format(suffix,fname)))

        bgSL = pyregion.parse("fk5; circle({0},{1},{2}\")"
                              .format(CL[0],
                                      CL[1],
                                      2*radius))
        bgsc = cube.subcube_from_ds9region(bgSL)
        print(bgsc)
        npix = np.count_nonzero(np.isfinite(bgsc[1000,:,:]))
        assert npix > 0
        bgsum = bgsc.sum(axis=(1,2))
        assert np.any(np.isfinite(bgsum))
        bgspec = (bgsum - spsum) / npix
        bgspec.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in scube.beams]),
                                              minor=np.nanmedian([bm.minor.to(u.deg).value for bm in scube.beams]),
                                              pa=np.nanmedian([bm.pa.to(u.deg).value for bm in scube.beams]),
                                             )
        bghdu = bgspec.hdu
        bghdu.header['OBJECT'] = name
        bghdu.writeto(spath("{0}_{1}_background_mean{2}.fits".format(fname,
                                                                     name,
                                                                     suffix)
                                ), clobber=True)

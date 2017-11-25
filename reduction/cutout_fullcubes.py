import spectral_cube
from astropy import units as u

for spw in range(4):
    cube = spectral_cube.SpectralCube.read('full_OrionSourceI_B3_spw{0}_lines.fits'.format(spw))

    sourceIsubcube = cube[:,2200:2600,2200:2600]

    sourceIsubcube.write('full_OrionSourceI_B3_spw{0}_lines_cutout.fits'.format(spw), overwrite=True)

    BNsubcube = cube[:,2850:2950,2730:2830]

    BNsubcube.write('full_OrionBN_B3_spw{0}_lines_cutout.fits'.format(spw), overwrite=True)

    if spw == 0:
        siov2cube = sourceIsubcube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=85.640452*u.GHz).spectral_slab(-40*u.km/u.s, 60*u.km/u.s)
        siov2cube.write('OrionSourceI_SiO_v=2_J=2-1_cutout.fits', overwrite=True)
        siov2cube.beam_threshold=0.1
        siov2cube.spectral_slab(-20*u.km/u.s, 35*u.km/u.s).moment1().write('OrionSourceI_SiO_v=2_J=2-1_moment1.fits', overwrite=True)

        sio29v0cube = sourceIsubcube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=85.759188*u.GHz).spectral_slab(-40*u.km/u.s, 60*u.km/u.s)
        sio29v0cube.write('OrionSourceI_29SiO_v=0_J=2-1_cutout.fits', overwrite=True)
        sio29v0cube.beam_threshold=0.1
        sio29v0cube.spectral_slab(-20*u.km/u.s, 35*u.km/u.s).moment1().write('OrionSourceI_29SiO_v=0_J=2-1_moment1.fits', overwrite=True)

        sio30v0cube = sourceIsubcube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.846995*u.GHz).spectral_slab(-40*u.km/u.s, 60*u.km/u.s)
        sio30v0cube.write('OrionSourceI_SiO_v=0_J=2-1_cutout.fits', overwrite=True)
        sio30v0cube.beam_threshold=0.1
        sio30v0cube.spectral_slab(-20*u.km/u.s, 35*u.km/u.s).moment1().write('OrionSourceI_30SiO_v=0_J=2-1_moment1.fits', overwrite=True)

        sio30v1cube = sourceIsubcube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.24343*u.GHz).spectral_slab(-40*u.km/u.s, 60*u.km/u.s)
        sio30v1cube.write('OrionSourceI_SiO_v=1_J=2-1_cutout.fits', overwrite=True)
        sio30v1cube.beam_threshold=0.1
        sio30v1cube.spectral_slab(-20*u.km/u.s, 35*u.km/u.s).moment1().write('OrionSourceI_30SiO_v=1_J=2-1_moment1.fits', overwrite=True)

for spw in range(4):
    cube = spectral_cube.SpectralCube.read('full_OrionSourceI_B6_spw{0}_lines.fits'.format(spw))

    sourceIsubcube = cube[:, 3582-625:3582+625, 3560-625:3560+625]

    sourceIsubcube.write('full_OrionSourceI_B6_spw{0}_lines_cutout.fits'.format(spw), overwrite=True)

    BNsubcube = cube[:, 5560-300:5560+300, 5097-300:5097+300]

    BNsubcube.write('full_OrionBN_B6_spw{0}_lines_cutout.fits'.format(spw), overwrite=True)

import lines
from astropy import units as u
from spectral_cube import SpectralCube

for cubefn in [
    'OrionSourceI_only.B6.robust0.5.longbaselines.spw1.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B6.robust0.5.longbaselines.spw3.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B6.robust0.5.longbaselines.spw0.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B6.robust0.5.longbaselines.spw2.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B3.robust0.5.spw3.clarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B3.robust0.5.spw1.clarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B3.robust0.5.spw2.clarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B3.robust0.5.spw0.clarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B7.lb.robust0.5.spw0.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B7.lb.robust0.5.spw3.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B7.lb.robust0.5.spw2.maskedclarkclean10000_medsub.image.pbcor.fits',
    'OrionSourceI_only.B7.lb.robust0.5.spw1.maskedclarkclean10000_medsub.image.pbcor.fits',
]:

    cube = SpectralCube.read(cubefn)

    for line, freq in lines.disk_lines.items():
        if 'Cl' in line:
            scube = (cube
                     .with_spectral_unit(u.km/u.s, rest_value=freq,
                                         velocity_convention='radio')
                     .spectral_slab(-20*u.km/u.s, 30*u.km/u.s))

            scube.write('{0}_{1}'.format(line, cubefn))

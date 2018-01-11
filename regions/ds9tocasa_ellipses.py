import regions
from astropy import units as u

def write_crtf(regions, filename):
    with open(filename, 'w') as fh:

        fh.write('#CRTF\n')

        for reg in regions:
            major,minor = reg.width, reg.height
            pa = reg.angle + 90*u.deg
            if major < minor:
                major,minor = minor,major
                pa += 90*u.deg
            fh.write("ellipse [[{0}deg, {1}deg], [{2}arcsec, {3}arcsec], {4}deg] coord={5}\n"
                     .format(reg.center.ra.deg[0],
                             reg.center.dec.deg[0],
                             major.to(u.arcsec).value,
                             minor.to(u.arcsec).value,
                             pa.to(u.deg).value,
                             reg.center.name.upper(),
                            ))

if __name__ == "__main__":

    regfile = regions.read_ds9('deepcleanregions.reg')

    write_crtf(regfile, 'deepcleanregions.crtf')

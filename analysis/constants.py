from astropy import units as u
from astropy import coordinates
import regions
import pvextractor
import paths
from line_point_offset import offset_to_point

d_orion = 415*u.pc

vcen = 5.5*u.km/u.s




# constants used for centroid_planes and seifried_analysis

source = 'sourceI'
#diskycoord_list = pyregion.open(paths.rpath("{0}_disk_pvextract.reg"
#                                            .format(source)))[0].coord_list
#diskycoords = coordinates.SkyCoord(["{0} {1}".format(diskycoord_list[jj],
#                                                     diskycoord_list[jj+1])
#                                    for jj in range(0,
#                                                    len(diskycoord_list),
#                                                    2)], unit=(u.deg,
#                                                               u.deg),
#                                   frame='fk5')
diskycoord_list = regions.read_ds9(paths.rpath("{0}_disk_pvextract.reg"
                                               .format(source)))
diskycoords = coordinates.SkyCoord([diskycoord_list[0].start,
                                    diskycoord_list[0].end])

#source = coordinates.SkyCoord(83.81048617*u.deg, -5.37516858*u.deg, frame='icrs')
source = coordinates.SkyCoord(regions.read_ds9(paths.rpath('sourceI_center.reg'))[0].center)

extraction_path = pvextractor.Path(diskycoords, width=0.01*u.arcsec)
origin = offset_to_point(source.ra.deg,
                         source.dec.deg,
                         extraction_path)*u.deg

# only approximate
central_freqs = {'B6': 224*u.GHz, 'B7': 345*u.GHz, 'B3': 95*u.GHz}

"""
http://docs.astropy.org/en/stable/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch
"""
from astropy import log
from astropy.io import fits
from astropy import wcs
import numpy as np
import glob
import re
import os
from astropy.utils.console import ProgressBar
from spectral_cube import SpectralCube

nchans_total = {0: 1920, 1: 1920, 2: 1920, 3: 1920}
min_nchans = 1900
frange = {0: sorted([344052.746e6, 344052.746e6+976.562e3*1920]),
          1: sorted([346052.686e6, 346052.686e6+976.562e3*1920]),
          2: sorted([335818.518e6, 335818.518e6-976.562e3*1920]),
          3: sorted([333927.125e6, 333927.125e6-976.562e3*1920]),
         }
fstep = {0:976623.21569824, # Hz
         1:976623.21569824, # Hz
         2:976623.21569824, # Hz
         3:976623.21569824, # Hz
        }

# Extract the appropriate pixel indices from the file name.
# A more sophisticated approach is probably better, in which the individual
# cubes are inspected for their start/end frequencies.
# But, on the other hand, for this process to make any sense at all, you
# have to have done the original cube imaging right
def getinds(fn):
    inds = re.search('lines([0-9]*)-([0-9]*)', fn).groups()
    return [int(ii) for ii in inds]

def get_max_ind(globstr):
    # replace nchans_total with the correct version from the actual data on
    # disk
    files = glob.glob(globstr)
    if len(files) == 0:
        return -1
    maxind = max([max(getinds(fn)) for fn in files])
    return maxind

def make_spw_cube(spw='spw{0}', spwnum=0, fntemplate='OrionSourceI',
                  overwrite_existing=False, bmaj_limits=None,
                  fnsuffix="", filesuffix='image.pbcor.fits',
                  first_endchannel='*',
                  cropends=False,
                  minimize=True,
                  debug_mode=False,
                  check_last_plane=False,
                  add_beam_info=True):
    """
    Parameters
    ----------
    spw : str
        String template for the input/output name
    spwnum : int
        The spectral window number
    fntemplate : str
        Filename template (goes into the glob)
    overwrite_existing : bool
        Overwrite data in the output cube?
    cropends: bool or int
        Number of pixels to crop off the ends of an image
    minimize: bool
        Compute the spatial minimal subcube before building the cube?  Slices
        for all subsequent cubes will be computed from the first cube.
    """
    if debug_mode:
        lvl = log.getEffectiveLevel()
        log.setLevel('DEBUG')

    spw = spw.format(spwnum)

    big_filename = '{1}_{0}{2}_lines.fits'.format(spw, fntemplate, fnsuffix)

    header_fn = glob.glob('OrionSourceI.B7.{0}.lines0-{4}.maskedclarkclean1000.{3}'
                          .format(spw, fntemplate, fnsuffix, filesuffix,
                                  first_endchannel))
    if len(header_fn) != 1:
        raise ValueError("Found too many or too few matches: {0}".format(header_fn))
    else:
        header_fn = header_fn[0]

    # First set up an empty file
    if not os.path.exists(big_filename):
        log.info("Creating large cube based on header {0}".format(header_fn))

        if minimize:
            cube0 = SpectralCube.read(header_fn)
            slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                    spatial_only=True)
            # use the calculated 3rd dimension, plus the difference of the
            # x and y slices
            #header['NAXIS2'] = slices[1].stop-slices[1].start
            #header['NAXIS1'] = slices[2].stop-slices[2].start
            header = cube0[slices].header
        else:
            header = fits.getheader(header_fn)

        # Make an arbitrary, small data before prepping the header
        data = np.zeros((100, 100), dtype=np.float32)
        hdu = fits.PrimaryHDU(data=data, header=header)
        cdelt_sign = np.sign(hdu.header['CDELT3'])
        # Set the appropriate output size (this can be extracted from the LISTOBS)
        naxis3_in = header['NAXIS3']
        header['NAXIS3'] = nchans_total[spwnum]
        header_wcs = wcs.WCS(fits.getheader(header_fn))
        header_specwcs = header_wcs.sub([wcs.WCSSUB_SPECTRAL])
        if cdelt_sign == -1:
            ind0, ind1 = getinds(header_fn)
            #5/20/2017: redoing some of this, and the text below is frightening but no longer relevant
            # a +1 was on the next line before an edit on 4/10/2017
            # it may have been rendered irrelevant when I included +1
            # channel in each cube?  Not clear - the arithmetic no longer
            # makes sense but is empirically necessary.
            assert ind0 == 0

            # these reindex the cube so that it has an increasing cdelt.
            header['CRPIX3'] = 1 #nchans_total[spwnum]
            header['CRVAL3'] = header_specwcs.wcs_pix2world([nchans_total[spwnum]],1)[0][0]
            header['CDELT3'] = np.abs(header_specwcs.wcs.cdelt[0])

            # ensure that the new CRVAL evaluated at its own position matches
            # the CRVAL3.  This should be impossible to fail unless WCS itself
            # fails
            newheaderspecwcs = wcs.WCS(header).sub([wcs.WCSSUB_SPECTRAL])
            crval3 = newheaderspecwcs.wcs_pix2world([header['CRPIX3']], 1)[0][0]
            np.testing.assert_array_almost_equal_nulp(crval3, header['CRVAL3'])


        shape = (header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])



        # Write to disk
        header.tofile(big_filename)
        # Using the 'append' io method, update the *header*
        with open(big_filename, 'rb+') as fobj:
            # Seek past the length of the header, plus the length of the
            # data we want to write.
            # The -1 is to account for the final byte that we are about to
            # write:
            # 'seek' works on bytes, so divide #bits / (bytes/bit)
            fobj.seek(len(header.tostring()) + (shape[0] *
                                                shape[1] *
                                                shape[2] *
                                                int(np.abs(header['BITPIX'])/8)) -
                      1)
            fobj.write(b'\0')

        big_cube = SpectralCube.read(big_filename)
        header_cube = SpectralCube.read(header_fn)
        # in both cases, SpectralCube sorts the extrema
        if cdelt_sign == 1:
            np.testing.assert_array_almost_equal_nulp(big_cube.spectral_extrema[0].value,
                                                      header_cube.spectral_extrema[0].value)
            np.testing.assert_array_almost_equal_nulp(big_cube.wcs.wcs.cdelt,
                                                      header_cube.wcs.wcs.cdelt)
        elif cdelt_sign == -1:
            np.testing.assert_array_almost_equal_nulp(big_cube.spectral_extrema[1].value,
                                                      header_cube.spectral_extrema[1].value)
            np.testing.assert_array_almost_equal_nulp(big_cube.wcs.wcs.cdelt[-1]*-1,
                                                      header_cube.wcs.wcs.cdelt[-1])

        log.info("Cube creation completed.  Now moving on to populating it.")


    # Find the appropriate files (this is NOT a good way to do this!  Better to
    # provide a list.  But wildcards are quick & easy...
    fileglob = "OrionSourceI.B7.{0}.lines*{3}".format(spw, fntemplate, fnsuffix,
                                                  filesuffix)
    files = glob.glob(fileglob)
    log.info("Files to be merged with glob {0}: ".format(fileglob))
    log.info(str(files))

    # open the file in update mode (it should have the right dims now)
    hdul = fits.open(big_filename, mode='update')
    main_wcs = wcs.WCS(hdul[0].header).sub([wcs.WCSSUB_SPECTRAL])

    if add_beam_info:
        shape = hdul[0].data.shape[0]
        if len(hdul) > 1 and isinstance(hdul[1], fits.BinTableHDU):
            pass
        else:
            hdul.append(fits.BinTableHDU(np.recarray(shape,
                                                     names=['BMAJ','BMIN','BPA','CHAN','POL'],
                                                     formats=['f4','f4','f4','i4','i4'])))

    # sorted so that we deal with zero first, since it has potential to be a problem.
    for fn in ProgressBar(sorted(files)):
        log.info("inds={0} fn={1}".format(getinds(fn), fn))
        ind0,ind1 = getinds(fn)

        # this is not correct...?
        # or maybe it only applies if cropends is set....
        # if ind0 == 0:
        #     ind1 = ind1 + 1

        cdelt = fits.getheader(fn)['CDELT3']
        if 'cdelt_sign' not in locals():
            cdelt_sign = np.sign(cdelt)
            log.warn("cdelt_sign was not defined: overwriting a"
                     " previously-existing file.  "
                     "This may not be what you want; the data could be going "
                     "opposite the parent cube.  Check that the original "
                     "header is OK. sign(CDELT) is now {0}, "
                     "while for the big header it is {1}"
                     .format(cdelt_sign,
                             np.sign(fits.getheader(big_filename)['CDELT3'])))

        if cropends:
            # don't crop 1st or last pixel in full cube
            if ind0 > 0:
                log.debug("ind0 going from {0} to {1}".format(ind0,ind0+cropends))
                ind0 = ind0 + cropends
                if cdelt_sign == 1:
                    dataind0 = cropends
                    log.debug("dataind0 going to {0}".format(cropends))
                else:
                    dataind1 = -cropends
                    log.debug("dataind1 going to {0}".format(-cropends))
            else:
                if cdelt_sign == 1:
                    dataind0 = 0
                    log.debug("dataind0 going to {0}".format(0))
                elif cdelt_sign == -1:
                    log.debug("dataind1 going to {0}".format(None))
                    dataind1 = None

            if (ind1 < nchans_total[spwnum] - 1):
                log.debug("ind1 going from {0} to {1}".format(ind1,ind1-cropends))
                ind1 = ind1 - cropends
                if cdelt_sign == 1:
                    dataind1 = - cropends
                    log.debug("dataind1 going to {0}".format(-cropends))
                elif cdelt_sign == -1:
                    dataind0 = cropends
                    log.debug("dataind0 going to {0}".format(cropends))
            else:
                if cdelt_sign == 1:
                    dataind1 = None
                else:
                    log.debug("dataind0 going to {0}".format(0))
                    dataind0 = 0
        else:
            dataind0 = 0
            dataind1 = None

        if cdelt_sign == -1:
            log.debug("Reversing indices from {0} {1} to ".format(ind0,ind1))
            ind1, ind0 = (nchans_total[spwnum] - ind0,
                          nchans_total[spwnum] - ind1)
            log.debug("{0} {1}".format(ind0, ind1))
            if ind0 < 0:
                ind0 = 0

        log.info("inds have been remapped to {0}, {1}".format(ind0, ind1))


        plane = hdul[0].data[ind0]
        lastplane = hdul[0].data[ind1-1]
        if np.all(plane == 0) or overwrite_existing or (check_last_plane and np.all(lastplane==0)):
            log.info("Replacing indices {0}->{2} {1}"
                     .format(getinds(fn), fn, (ind0,ind1)))

            data = fits.getdata(fn)
            dwcs = wcs.WCS(fits.getheader(fn)).sub([wcs.WCSSUB_SPECTRAL])

            dataind1 = data.shape[0]+(dataind1 or 0)

            # handle the case where I made the indices NOT match the cube...
            # this is really stupid and should be removed because I should have
            # made the input cubes correct.  Oh well.
            if np.abs(ind1 - ind0) < np.abs(dataind1 - dataind0):
                dataind1 = dataind0 + np.abs(ind1-ind0)

            if cdelt_sign == -1:
                dataind0, dataind1 = dataind1, dataind0
                dwcs0 = dwcs.wcs_pix2world([dataind0-1], 0)[0][0]
                dwcs1 = dwcs.wcs_pix2world([dataind1], 0)[0][0]
            else:
                dwcs0 = dwcs.wcs_pix2world([dataind0], 0)[0][0]
                dwcs1 = dwcs.wcs_pix2world([dataind1-1], 0)[0][0]
            hwcs0 = main_wcs.wcs_pix2world([ind0], 0)[0][0]
            hwcs1 = main_wcs.wcs_pix2world([ind1-1], 0)[0][0]
            
            if not np.isclose(hwcs0, dwcs0, atol=0.5*np.abs(cdelt), rtol=0):
                log.error("current data, big cube indices: {0},{1} and {2},{3}"
                          .format(dataind0,dataind1,ind0,ind1))
                raise ValueError("World coordinates of first pixels do not match: {0} - {1} = {2} ({3} cdelt)"
                                 .format(dwcs0,hwcs0,dwcs0-hwcs0,(dwcs0-hwcs0)/cdelt))
            if not np.isclose(hwcs1, dwcs1, atol=0.5*np.abs(cdelt), rtol=0):
                log.error("current data, big cube indices: {0},{1} and {2},{3}"
                          .format(dataind0,dataind1,ind0,ind1))
                raise ValueError("World coordinates of last pixels do not match: {0} - {1} = {2} ({3} cdelt)"
                                 .format(dwcs1,hwcs1,dwcs1-hwcs1,(dwcs1-hwcs1)/cdelt))

            if 'slices' not in locals():
                if minimize:
                    log.info("Determining slices")
                    cube0 = SpectralCube.read(header_fn)
                    slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                            spatial_only=True)
                    log.info("Slices are {0}".format(slices))
                else:
                    slices = (slice(None),)*3


            if bmaj_limits is not None:
                log.info("Identifying acceptable beams")
                beamtable = fits.open(fn)[1]
                ok_beam = ((beamtable.data['BMAJ'] > bmaj_limits[0]) &
                           (beamtable.data['BMAJ'] < bmaj_limits[1]))
                data[~ok_beam] = np.nan
                log.info("Found {0} bad beams of {1}".format((~ok_beam).sum(),
                                                             ok_beam.size))

            if cdelt_sign == -1:
                if dataind1 == 0:
                    dataslice = slice(dataind0-1, None, -1)
                elif dataind1 >= 1:
                    dataslice = slice(dataind0-1, dataind1-1, -1)
                else:
                    raise ValueError("Something is wrong with dataind0")
            else:
                dataslice = slice(dataind0, dataind1, 1)
            log.info("Dataslice is {0}".format(dataslice))

            assert hdul[0].data[ind0:ind1].shape == data[dataslice, slices[1], slices[2]].shape

            if not debug_mode:
                if add_beam_info:
                    log.info("Adding beam information")
                    beamtable = fits.open(fn)[1]
                    hdul[1].data[ind0:ind1] = beamtable.data[dataslice]


                log.info("Inserting data")
                hdul[0].data[ind0:ind1,:,:] = data[dataslice, slices[1], slices[2]]
                log.info("Flushing")
                hdul.flush()
                log.info("Done with iteration for {0}".format(fn))

    if debug_mode:
        log.setLevel(lvl)

if __name__ == "__main__":
    for spw in (0,1,2,3):

        mxind = get_max_ind('OrionSourceI.B7.spw{0}.lines*fits'.format(spw, ))
        if mxind < min_nchans:
            log.critical("Skipping {0} b/c only {1} chans".format(spw, mxind))
            continue
        nchans_total[spw] = mxind
        log.info("nchans_total[{0}] = {1}".format(spw, mxind))

        if os.path.exists('OrionSourceI.B7.spw{0}.lines0-60.maskedclarkclean1000.image.pbcor.fits'.format(spw)):

            make_spw_cube(spw='spw{0}', spwnum=spw,
                          fntemplate='full_OrionSourceI_B7_maskedclean',
                          overwrite_existing=False, bmaj_limits=None,
                          fnsuffix="", filesuffix='image.pbcor.fits',
                          first_endchannel=60,
                          check_last_plane=True,
                          #debug_mode=True,
                          cropends=0, minimize=True, add_beam_info=True)

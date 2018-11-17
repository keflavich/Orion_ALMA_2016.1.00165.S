import numpy as np
import os
from astropy.io import fits
from astropy import wcs
import regions

from taskinit import iatool, casalog

from importfits_cli import importfits_cli as importfits
from exportfits_cli import exportfits_cli as exportfits
from makemask_cli import makemask_cli as makemask

ia = iatool()

def get_mask(regions, hdu):

    ww = wcs.WCS(hdu.header).celestial

    reg_pix = [reg.to_pixel(ww) for reg in regions]

    reg_masks = [reg.to_mask() for reg in reg_pix]

    mask = np.zeros(hdu.data.shape[-2:], dtype='int16')

    assert mask.size > 1

    for mm,reg in zip(reg_masks, regions):
        if mm.cutout(mask) is None:
            #casalog.post("Skipped region {0}".format(str(reg)), origin='get_mask')
            continue
        else:
            #nmask = mask.sum()
            if hasattr(mm, '_to_image_partial_overlap'):
                mm._to_image_partial_overlap(mask)
            elif hasattr(mm, 'to_image_partial_overlap'):
                mm.to_image_partial_overlap(mask)
            else:
                raise ValueError("Something is wrong with the regions build")
            #nmask_after = mask.sum()
            #casalog.post("For region {0}, mask went from {1} to {2}"
            #             .format(str(reg), nmask, nmask_after),
            #             origin='get_mask',
            #            )

    nmask_after = mask.sum()
    casalog.post("Final mask in get_mask has {0} of {1} pixels included ({2}%)"
                 .format(nmask_after, mask.size, nmask_after/float(mask.size)*100),
                 origin='get_mask',
                )

    return mask

def reg_to_mask(regfile, baseimage):
    """
    This is like makemask(mode='copy', inpimage=baseimage, inpmask=regfile)

    """
    casalog.post("Creating mask from regions {0} into base image {1}".format(regfile, baseimage), origin='reg_to_mask')

    assert baseimage.endswith('.image') or baseimage.endswith('.image.tt0')
    cleanbox_mask = baseimage.replace(".image.tt0", ".image").replace(".image", ".mask")
    cleanbox_mask_image = baseimage

    if not os.path.exists(cleanbox_mask) or not os.path.exists(cleanbox_mask_image):

        # create a mask based on region selection (no thresholding here)
        dirtyimagename = baseimage

        exportfits(dirtyimagename, dirtyimagename+".fits", overwrite=True)
        reg = regions.read_ds9(regfile)
        imghdu = fits.open(dirtyimagename+".fits")[0]
        assert imghdu.data.shape[-1] != 1
        assert imghdu.data.shape[-2] != 1

        mask = get_mask(reg, imghdu)
        assert mask.shape[0] > 1
        assert mask.shape[1] > 1
        assert mask.size > 1

        imghdu.data = mask.astype('int16')
        imghdu.header['BITPIX'] = 16
        imghdu.writeto(cleanbox_mask+'.fits', clobber=True)
        importfits(fitsimage=cleanbox_mask+'.fits',
                   imagename=cleanbox_mask_image,
                   overwrite=True)
        ia.open(cleanbox_mask_image)
        ia.calcmask(mask="'{0}' > 0.5".format(cleanbox_mask_image),
                    name='cleanbox_mask')

        ia.close()
        makemask(mode='copy', inpimage=cleanbox_mask_image,
                 inpmask=cleanbox_mask_image+":cleanbox_mask",
                 output=cleanbox_mask,
                 overwrite=True)

    mask = cleanbox_mask

    ia.open(mask)
    stats = ia.statistics()
    ia.close()

    casalog.post("Resulting mask file is {0}".format(mask), origin='reg_to_mask')
    casalog.post("Mask sum = {0} out of {1}, or {2}%".format(stats['sum'],
                                                             stats['npts'],
                                                             stats['sum']/stats['npts']*100),
                 origin='reg_to_mask')

    return mask

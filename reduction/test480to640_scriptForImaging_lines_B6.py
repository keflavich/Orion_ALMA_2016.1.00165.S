
import numpy as np
import os
import glob
import datetime

def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.psf', fitsimage=myimagebase+'.psf.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.sumwt', fitsimage=myimagebase+'.sumwt.fits', dropdeg=True, overwrite=True)



mslist = ['uid___A001_X88e_X1d3_calibrated.ms',
          'uid___A002_Xb925ef_X4334_calibrated.ms']

for spw,spws in enumerate([(0,), ]):

    for suffix, niter in (('clarkclean1000', 1000), ('dirty', 0)):
        
        step = 1920/32
        for startchan in np.arange(0, 1920, step):

            if startchan not in (480, 540):
                continue

            imagename = 'OrionSourceI.B6.spw{0}.lines{2}-{3}.{1}'.format(spw, suffix, startchan, startchan+step)

            # DO NOT skip any, just redo
            #if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
            #    print("Skipping completed file {0}".format(imagename))
            #    continue

            print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
            tclean(vis=mslist,
                   imagename=imagename,
                   datacolumn='data',
                   spw=",".join(['{0}'.format(ss) for ss in spws]),
                   field='Orion_BNKL_source_I',
                   specmode='cube',
                   outframe='LSRK',
                   start=startchan,
                   nchan=step,
                   threshold='1mJy',
                   imsize=[7168, 7168],
                   cell=['0.004arcsec'],
                   niter=niter,
                   deconvolver='clark',
                   gridder='standard',
                   weighting='briggs',
                   robust=0.5,
                   pbcor=True,
                   pblimit=0.2,
                   savemodel='none',
                   parallel=True,
                   interactive=False)
            makefits(imagename)

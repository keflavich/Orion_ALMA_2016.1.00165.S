from astropy.io import fits
import numpy as np
import os
import glob
import datetime

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

    if cleanup:
        for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                       'mask', 'image', 'workdirectory'):
            os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))




mslist = ['member.uid___A001_X88e_X1d9_calibrated.ms']

for ms in mslist:
    listobs(ms, listfile=ms+'.listobs', overwrite=True)

for spw,spws in enumerate([(0,4), (1,5), (2,6), (3,7)]):

    for suffix, niter in (('clarkclean10000', 10000), ):

        for robust in (-2, 0.5, 2):

            imagename = 'OrionSourceI_only.B3.robust{2}.spw{0}.{1}'.format(spw, suffix, robust)


            if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
                    print("Skipping completed file {0}".format(imagename))
                    continue
                else:
                    print("Redoing {0} because it's in fk5.".format(imagename))

            print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
            tclean(vis=mslist,
                   imagename=imagename,
                   datacolumn='data',
                   spw=",".join(['{0}'.format(ss) for ss in spws]),
                   field='Orion_BNKL_source_I',
                   specmode='cube',
                   outframe='LSRK',
                   threshold='15mJy',
                   imsize=[128, 128],
                   cell=['0.008arcsec'],
                   niter=niter,
                   cycleniter=-1, # -1 is default
                   cyclefactor=0.0001, # set very small: try to prevent major cycles
                   phasecenter='ICRS 5h35m14.5184s -5d22m30.6199s',
                   deconvolver='clark',
                   gridder='standard',
                   weighting='briggs',
                   robust=robust,
                   pbcor=True,
                   pblimit=0.2,
                   savemodel='none',
                   chanchunks=1,
                   parallel=True,
                   interactive=False)
            makefits(imagename)


            imagename = 'OrionBN_only.B3.robust{2}.spw{0}.{1}'.format(spw, suffix, robust)
            if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                if fits.getheader("{0}.image.pbcor.fits".format(imagename))['RADESYS'] == 'ICRS':
                    print("Skipping completed file {0}".format(imagename))
                    continue
                else:
                    print("Redoing {0} because it's in fk5.".format(imagename))

            print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
            tclean(vis=mslist,
                   imagename=imagename,
                   datacolumn='data',
                   spw=",".join(['{0}'.format(ss) for ss in spws]),
                   field='Orion_BNKL_source_I',
                   specmode='cube',
                   outframe='LSRK',
                   threshold='15mJy',
                   imsize=[128, 128],
                   cell=['0.008arcsec'],
                   niter=niter,
                   cycleniter=-1, # -1 is default
                   cyclefactor=0.0001, # set very small: try to prevent major cycles
                   phasecenter='ICRS 5h35m14.108 -5d22m22.669s',
                   deconvolver='clark',
                   gridder='standard',
                   weighting='briggs',
                   robust=robust,
                   pbcor=True,
                   pblimit=0.2,
                   savemodel='none',
                   chanchunks=1,
                   parallel=True,
                   interactive=False)
            makefits(imagename)

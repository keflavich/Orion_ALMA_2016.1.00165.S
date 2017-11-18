
import numpy as np
import os
import glob
import datetime

def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image



mslist = ['member.uid___A001_X88e_X1d9_calibrated.ms']

for ms in mslist:
    listobs(ms, listfile=ms+'.listobs', overwrite=True)

for spw,spws in enumerate([(0,4), (1,5), (2,6), (3,7)]):

    for suffix, niter in (('clarkclean1000', 1000), ):

        step = 1920/32
        for startchan in np.arange(0, 1920, step):

            imagename = 'OrionSourceI.B3.spw{0}.lines{2}-{3}.{1}'.format(spw, suffix, startchan, startchan+step)
            if not os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
                tclean(vis=mslist,
                       imagename=imagename,
                       datacolumn='data',
                       spw=",".join(['{0}'.format(ss) for ss in spws]),
                       field='Orion_BNKL_source_I',
                       specmode='cube',
                       outframe='LSRK',
                       threshold='1mJy',
                       imsize=[4800, 4800],
                       cell=['0.016arcsec'],
                       niter=niter,
                       deconvolver='clark',
                       gridder='standard',
                       weighting='briggs',
                       robust=0.5,
                       pbcor=True,
                       pblimit=0.2,
                       savemodel='none',
                       chanchunks=1,
                       start=startchan,
                       nchan=step,
                       parallel=True,
                       interactive=False)
                makefits(imagename)

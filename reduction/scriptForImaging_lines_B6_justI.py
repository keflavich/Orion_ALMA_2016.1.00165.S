from astropy.io import fits

import numpy as np
import os
import glob
import datetime

from source_ids import sources_fmtd

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.psf', fitsimage=myimagebase+'.psf.fits', dropdeg=True, overwrite=True) # export the PSF image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

    if cleanup:
        for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                       'mask', 'image', 'workdirectory'):
            os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))



mslist = ['uid___A001_X88e_X1d3_calibrated.ms',
          'uid___A002_Xb925ef_X4334_calibrated.ms']

for ms in mslist:
    print("Listing {0}".format(ms))
    result = listobs(ms, listfile=ms+'.listobs', overwrite=True)
    print("Done listing {0}".format(ms))
    assert result,"Listing {0} failed".format(ms)

for sourcename, coordinate in sources_fmtd.items():

    for spw,spws in enumerate([(0,), (1,), (2,), (3,)]):

        for suffix, niter in (('maskedclarkclean10000', 10000), ):
            for robust in (-2, 0.5, 2):

                imagename = 'Orion{3}_only.B6.robust{2}.spw{0}.{1}'.format(spw, suffix, robust, sourcename)

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
                       cycleniter=-1, # -1 is default
                       cyclefactor=0.0001, # set very small: try to prevent major cycles
                       # clean only Source I and BN
                       #mask=['circle[[5h35m14.5184s,-5d22m30.6199s],0.222arcsec]'],
                       phasecenter=coordinate,
                       imsize=[128, 128],
                       cell=['0.004arcsec'],
                       niter=niter,
                       deconvolver='clark',
                       gridder='standard',
                       weighting='briggs',
                       robust=robust,
                       uvrange='10~5000000m',
                       pbcor=True,
                       pblimit=0.2,
                       savemodel='none',
                       parallel=True,
                       interactive=False)
                makefits(imagename)


            for robust in (-2, 0.5, 2):
                imagename = 'Orion{3}_only.B6.robust{2}.longbaselines.spw{0}.{1}'.format(spw, suffix, robust, sourcename)

                if os.path.exists("{0}.image.pbcor.fits".format(imagename)):
                    print("Skipping completed file {0}".format(imagename))
                    continue

                print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
                tclean(vis=mslist,
                       imagename=imagename,
                       datacolumn='data',
                       spw=",".join(['{0}'.format(ss) for ss in spws]),
                       field='Orion_BNKL_source_I',
                       specmode='cube',
                       outframe='LSRK',
                       threshold='15mJy',
                       cycleniter=-1, # -1 is default
                       cyclefactor=0.0001, # set very small: try to prevent major cycles
                       # clean only Source I and BN
                       #mask=['circle[[5h35m14.5184s,-5d22m30.6199s],0.222arcsec]'],
                       phasecenter=coordinate,
                       imsize=[128, 128],
                       cell=['0.004arcsec'],
                       uvrange='50~5000000m',
                       niter=niter,
                       deconvolver='clark',
                       gridder='standard',
                       weighting='briggs',
                       robust=robust,
                       pbcor=True,
                       pblimit=0.2,
                       savemodel='none',
                       parallel=True,
                       interactive=False)
                makefits(imagename)

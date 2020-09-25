import sys
sys.path.insert(0,'.')
from makemask_regions import reg_to_mask
import os

parallel = False

redo = False

def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt0', fitsimage=myimagebase+'.image.tt0.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.psf.tt0', fitsimage=myimagebase+'.psf.tt0.fits', dropdeg=True, overwrite=True) # export the PSF image


cell='0.004arcsec' # cell size for imaging.
imsize = [14000,14000] # size of image in pixels.
imsize = [30000,30000] # 14k gets you to 0.6 pb!  

#contvis = ['band6_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
#if not os.path.exists(contvis[0]):
#    make_cont_ms()
#contvis = 'B6_calibrated_final_cont.ms'

contvis = 'B3_calibrated_final_cont.ms'

for uvrange, uvrangename in (
                             ('10~300000m', 'allbaselines.huge'),
):
    for robust, depth1, depth2 in (
                                   (0.5, '1mJy', '0.05mJy',),
                                  ):

        contimagename = 'Orion_SourceI_B3_continuum_r{0}.clean{1}.{2}'.format(robust, depth1, uvrangename)
        if redo or not os.path.exists(contimagename+".residual.tt0"):
            # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
            os.system('rm -rf ' + contimagename + "*")
            tclean(vis=contvis,
                   imagename=contimagename,
                   datacolumn='data',
                   field='Orion_BNKL_source_I',
                   specmode='mfs',
                   deconvolver='mtmfs',
                   nterms=2,
                   scales=[0,4,12],
                   smallscalebias=0.8,
                   imsize = imsize,
                   cell= cell,
                   weighting = 'briggs',
                   robust = robust,
                   niter = int(1e6),
                   threshold = depth1,
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange=uvrange,
                   parallel=parallel,
                   pblimit=0.0,
                  )
            makefits(contimagename)

        prevcontimage = contimagename

        maskfile = reg_to_mask('deepcleanregions.reg', prevcontimage+".image.tt0")

        contimagename = 'Orion_SourceI_B3_continuum_r{0}.clean{1}.{2}.deepmask'.format(robust, depth2, uvrangename)
        if redo or not os.path.exists(contimagename+".image.tt0"):
            #os.system('rm -rf ' + contimagename + "*")
            if not os.path.exists(contimagename+".model.tt0"):
                startmodel=[prevcontimage+'.model.tt0', prevcontimage+'.model.tt1']
            else:
                startmodel=''
            tclean(vis=contvis,
                   imagename=contimagename,
                   datacolumn='data',
                   startmodel=startmodel,
                   field='Orion_BNKL_source_I',
                   specmode='mfs',
                   deconvolver='mtmfs',
                   mask=prevcontimage+".mask",
                   nterms=2,
                   scales=[0,4,12],
                   smallscalebias=0.8,
                   imsize = imsize,
                   cell= cell,
                   weighting = 'briggs',
                   robust = robust,
                   niter = int(1e6),
                   threshold = depth2,
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange=uvrange,
                   parallel=parallel,
                   pblimit=0.0,
                  )
            makefits(contimagename)

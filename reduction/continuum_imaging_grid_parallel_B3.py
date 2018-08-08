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

#contvis = ['band6_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
#if not os.path.exists(contvis[0]):
#    make_cont_ms()
#contvis = 'B6_calibrated_final_cont.ms'

contvis = 'B3_calibrated_final_cont.ms'

for uvrange, uvrangename in (('50~300000m', '50mplus'),
                             ('825~10000klambda', '825to10000kl'),
                             ('150~300000m', '150mplus'),
                             ('200~300000m', '200mplus'),
                             ('500~300000klambda', '500klplus'),
                             ('1000~300000klambda', '1000klplus'),
                             ('1500~300000klambda', '1500klplus'),
                             ('10~300000m', 'allbaselines'),):
    for robust, depth1, depth2 in ((2, '2mJy', '1mJy',),
                                   (0.5, '1mJy', '0.5mJy',),
                                   (-2, '0.5mJy', '0.1mJy')):

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
                   niter = int(1e5),
                   threshold = depth1,
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange=uvrange,
                   parallel=True,
                  )
            makefits(contimagename)

            prevcontimage = contimagename
            contimagename = 'Orion_SourceI_B3_continuum_r{0}.clean{1}.{2}.deepmask'.format(robust, depth2, uvrangename)
            os.system('rm -rf ' + contimagename + "*")
            tclean(vis=contvis,
                   imagename=contimagename,
                   datacolumn='data',
                   startmodel=[prevcontimage+'.model.tt0', prevcontimage+'.model.tt1'],
                   field='Orion_BNKL_source_I',
                   specmode='mfs',
                   deconvolver='mtmfs',
                   mask=['deepcleanregions.crtf'],
                   nterms=2,
                   scales=[0,4,12],
                   smallscalebias=0.8,
                   imsize = imsize,
                   cell= cell,
                   weighting = 'briggs',
                   robust = robust,
                   niter = int(1e5),
                   threshold = depth2,
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange=uvrange,
                   parallel=True,
                  )
            makefits(contimagename)

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
imsize = [7168,7168] # size of image in pixels.

#contvis = ['band6_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
#if not os.path.exists(contvis[0]):
#    make_cont_ms()
#contvis = 'B6_calibrated_final_cont.ms'

selfcal_vis = 'B6_selfcal.ms'
nocal_vis = 'B6_nocal.ms'

if not os.path.exists(nocal_vis):
    split(selfcal_vis, nocal_vis, datacolumn='data')

    clearcal(vis=nocal_vis)


for robust, depth1, depth2 in ((2, '2mJy', '1mJy',),
                               (0.5, '1mJy', '0.5mJy',),
                               (-2, '0.5mJy', '0.1mJy')):

    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.50mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=nocal_vis,
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
               uvrange='50~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.50mplus.deepmask'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=nocal_vis,
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
               uvrange='50~36000m',
               parallel=True,
              )
        makefits(contimagename)



    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.allbaselines'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=nocal_vis,
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
               uvrange='10~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.deepmask.allbaselines'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=nocal_vis,
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
               uvrange='10~36000m',
               parallel=True,
              )
        makefits(contimagename)





    applycal(vis=selfcal_vis, gaintable=["phase_4.cal"],
             interp="linear", applymode='calonly', calwt=False)


    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.50mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='50~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.deepmask.50mplus'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='50~36000m',
               parallel=True,
              )
        makefits(contimagename)



    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.allbaselines'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='10~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.deepmask.allbaselines'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='10~36000m',
               parallel=True,
              )
        makefits(contimagename)


    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.150mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
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
               uvrange='150~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.150mplus.deepmask'.format(robust, depth2)
        tclean(vis=selfcal_vis,
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
               uvrange='150~36000m',
               parallel=True,
              )
        makefits(contimagename)


    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.150mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='150~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.deepmask.150mplus'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='150~36000m',
               parallel=True,
              )
        makefits(contimagename)

    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.200mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='200~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.selfcal.phase4.deepmask.200mplus'.format(robust, depth2)
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
               imagename=contimagename,
               datacolumn='corrected',
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
               uvrange='200~36000m',
               parallel=True,
              )
        makefits(contimagename)


    contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.200mplus'.format(robust, depth1)
    if redo or not os.path.exists(contimagename+".residual.tt0"):
        # First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
        os.system('rm -rf ' + contimagename + "*")
        tclean(vis=selfcal_vis,
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
               uvrange='200~36000m',
               parallel=True,
              )
        makefits(contimagename)

        prevcontimage = contimagename
        contimagename = 'Orion_SourceI_B6_continuum_r{0}.clean{1}.200mplus.deepmask'.format(robust, depth2)
        tclean(vis=selfcal_vis,
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
               uvrange='200~36000m',
               parallel=True,
              )
        makefits(contimagename)



# what if we use the high-resolution data model as the input before cleaning the lower-resolution data?
prevcontimage = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.150mplus'
contimagename = 'Orion_SourceI_B6_continuum_r2.clean0.5mJy.selfcal.phase4.allbaselines.highresstartmodel'
if redo or not os.path.exists(contimagename+".residual.tt0"):
    os.system('rm -rf ' + contimagename + "*")
    tclean(vis=selfcal_vis,
           imagename=contimagename,
           datacolumn='corrected',
           startmodel=[prevcontimage+'.model.tt0', prevcontimage+'.model.tt1'],
           field='Orion_BNKL_source_I',
           specmode='mfs',
           deconvolver='mtmfs',
           nterms=2,
           scales=[0,4,12,36],
           imsize = imsize,
           cell= cell,
           weighting = 'briggs',
           robust = 2,
           niter = int(1e5),
           threshold = '2.5mJy',
           interactive = False,
           outframe='LSRK',
           veltype='radio',
           savemodel='none',
           uvrange='10~36000m',
           parallel=True,
          )
    makefits(contimagename)


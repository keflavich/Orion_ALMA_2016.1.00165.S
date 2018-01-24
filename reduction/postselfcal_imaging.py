raise ValueError("Use parallel version instead.")

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

clearcal(vis=selfcal_vis)



contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.50mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='50~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.50mplus.deepmask'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='50~36000m',
      )
makefits(contimagename)



contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.allbaselines'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='10~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.deepmask.allbaselines'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='10~36000m',
      )
makefits(contimagename)





applycal(vis=selfcal_vis, gaintable=["phase_4.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4.50mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='50~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.50mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='50~36000m',
      )
makefits(contimagename)



contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4.allbaselines'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='10~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='10~36000m',
      )
makefits(contimagename)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.150mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='150~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.150mplus.deepmask'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='150~36000m',
      )
makefits(contimagename)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4.150mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='150~36000m',
      )
makefits(contimagename)

prevcontimage = contimagename
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.150mplus'
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
       robust = -2,
       niter = int(1e5),
       threshold = '0.1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='none',
       uvrange='150~36000m',
      )
makefits(contimagename)

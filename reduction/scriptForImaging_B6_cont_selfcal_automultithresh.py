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
contvis = 'B6_calibrated_final_cont.ms'

selfcal_vis = 'B6_selfcal_automultithresh.ms'
os.system('rm -rf {0}'.format(selfcal_vis))
os.system('cp -r {0} {1}'.format(contvis, selfcal_vis))

flagdata(vis=selfcal_vis, mode='manual', autocorr=True, uvrange='0~12m')


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean5mJy.automultithresh'

for suffix in ('', '.tt0', '.tt1', '.tt2'):
    for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
        todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
        if os.path.exists(todel):
            os.system('rm -rf {0}'.format(todel))


tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12,48],
       smallscalebias=0.8,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       usemask='auto-multithresh',
       threshold = '5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )

makefits(contimagename)


rmtables(['automultithresh_phase_0.cal'])
gaincal(vis=selfcal_vis, caltable='automultithresh_phase_0.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('automultithresh_phase_0.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='automultithresh_phase_0_vs_time.png')

applycal(vis=selfcal_vis, gaintable=["automultithresh_phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean4mJy.automultithresh.selfcal.phase0'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12,48],
       smallscalebias=0.8,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       usemask='auto-multithresh',
       threshold = '4mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


rmtables(['automultithresh_phase_1.cal'])
gaincal(vis=selfcal_vis, caltable='automultithresh_phase_1.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('automultithresh_phase_1.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='automultithresh_phase_1_vs_time.png')

applycal(vis=selfcal_vis, gaintable=["automultithresh_phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean3mJy.automultithresh.selfcal.phase1'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12,48],
       smallscalebias=0.8,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       usemask='auto-multithresh',
       threshold = '3mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)



rmtables(['automultithresh_phase_2.cal'])
gaincal(vis=selfcal_vis, caltable='automultithresh_phase_2.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('automultithresh_phase_2.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='automultithresh_phase_2_vs_time.png')

applycal(vis=selfcal_vis, gaintable=["automultithresh_phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean2mJy.automultithresh.selfcal.phase2'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12,48],
       smallscalebias=0.8,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       usemask='auto-multithresh',
       threshold = '2mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


rmtables(['automultithresh_phase_3.cal'])
gaincal(vis=selfcal_vis, caltable='automultithresh_phase_3.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('automultithresh_phase_3.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='automultithresh_phase_3_vs_time.png')

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

selfcal_vis = 'B6_selfcal.ms'
os.system('rm -rf {0}'.format(selfcal_vis))
os.system('cp -r {0} {1}'.format(contvis, selfcal_vis))

flagdata(vis=selfcal_vis, mode='manual', autocorr=True, uvrange='0~12m')


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean5mJy'

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
       scales=[0,4,12],
       smallscalebias=0.8,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       threshold = '5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )

makefits(contimagename)


rmtables(['phase_0.cal'])
gaincal(vis=selfcal_vis, caltable='phase_0.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_0.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_0_vs_time.png')

clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean4mJy.selfcal.phase0'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       threshold = '4mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


rmtables(['phase_1.cal'])
gaincal(vis=selfcal_vis, caltable='phase_1.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_1.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_1_vs_time.png')

clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean3mJy.selfcal.phase1'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       threshold = '3mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)



rmtables(['phase_2.cal'])
gaincal(vis=selfcal_vis, caltable='phase_2.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_2.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_2_vs_time.png')

clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean2mJy.selfcal.phase2'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       threshold = '2mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


rmtables(['phase_3.cal'])
gaincal(vis=selfcal_vis, caltable='phase_3.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_3.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_3_vs_time.png')

clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)

contimagename = 'Orion_SourceI_B6_continuum_r-2.clean1mJy.selfcal.phase3'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       threshold = '1mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


rmtables(['phase_4.cal'])
gaincal(vis=selfcal_vis, caltable='phase_4.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_4.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_4_vs_time.png')

# compute, but don't use, amp self-cal
rmtables(['amp_4.cal'])
gaincal(vis=selfcal_vis, caltable='amp_4.cal', solint='30s', gaintype='G',
        calmode='a')

plotcal('amp_4.cal', xaxis='time', yaxis='amp', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='amp_4_vs_time.png')


clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_4.cal"],
         interp="linear", applymode='calonly', calwt=False)



contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask.allbaselines'
# First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       savemodel='modelcolumn',
       uvrange='10~36000m',
      )

tclean(vis=selfcal_vis,
       imagename=contimagename,
       startmodel=[contimagename+'.model.tt0', contimagename+'.model.tt1'],
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
       savemodel='modelcolumn',
       uvrange='10~36000m',
      )
makefits(contimagename)




contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.phase4.deepmask'
# First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )

tclean(vis=selfcal_vis,
       imagename=contimagename,
       startmodel=[contimagename+'.model.tt0', contimagename+'.model.tt1'],
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)




rmtables(['phase_5.cal'])
gaincal(vis=selfcal_vis, caltable='phase_5.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_5.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_5_vs_time.png')

rmtables(['amp_5.cal'])
gaincal(vis=selfcal_vis, caltable='amp_5.cal', solint='30s', gaintype='G',
        calmode='a')

plotcal('amp_5.cal', xaxis='time', yaxis='amp', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='amp_5_vs_time.png')

rmtables(['amp_5_longbaselines.cal'])
gaincal(vis=selfcal_vis, caltable='amp_5_longbaselines.cal', solint='30s', gaintype='G',
        uvrange='150~10000000m',
        calmode='a')

plotcal('amp_5_longbaselines.cal', xaxis='time', yaxis='amp', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='amp_5_longbaselines_vs_time.png')
plotcal('amp_5_longbaselines.cal', xaxis='antenna', yaxis='amp',
        figfile='amp_5_longbaselines_vs_antenna.png')

clearcal(vis=selfcal_vis)

applycal(vis=selfcal_vis, gaintable=["phase_5.cal", "amp_5_longbaselines.cal"],
         interp="linear", applymode='calonly', calwt=False)

contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)

# imaging experiments:
contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.ampphase5.deepmask'
# First, clean everything to 0.5 mJy/beam, then clean just the specified regions deeper
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )

tclean(vis=selfcal_vis,
       imagename=contimagename,
       startmodel=[contimagename+'.model.tt0', contimagename+'.model.tt1'],
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)


contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5.largescalefirst'
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
       threshold = '2.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
tclean(vis=selfcal_vis,
       imagename=contimagename,
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)

contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5.centeredBN'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       phasecenter='ICRS 5h35m14.108 -5d22m22.669s',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12],
       smallscalebias=0.8,
       imsize = [3000,3000],
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)




contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5.bias0.95'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       nterms=2,
       scales=[0,4,12,48],
       smallscalebias=0.95,
       imsize = imsize,
       cell= cell,
       weighting = 'briggs',
       robust = -2,
       niter = int(1e5),
       threshold = '0.5mJy',
       interactive = False,
       outframe='LSRK',
       veltype='radio',
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)





contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.ampphase5.wterms'
os.system('rm -rf ' + contimagename + "*")
tclean(vis=selfcal_vis,
       imagename=contimagename,
       field='Orion_BNKL_source_I',
       specmode='mfs',
       deconvolver='mtmfs',
       gridder='wproject',
       wprojplanes=30,
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
       savemodel='modelcolumn',
       uvrange='50~36000m',
      )
makefits(contimagename)





rmtables(['phase_6.cal'])
gaincal(vis=selfcal_vis, caltable='phase_6.cal', solint='int', gaintype='G',
        calmode='p')

plotcal('phase_6.cal', xaxis='time', yaxis='phase', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='phase_6_vs_time.png')

rmtables(['amp_6.cal'])
gaincal(vis=selfcal_vis, caltable='amp_6.cal', solint='30s', gaintype='G',
        calmode='a')

plotcal('amp_6.cal', xaxis='time', yaxis='amp', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='amp_6_vs_time.png')

rmtables(['amp_6_longbaselines.cal'])
gaincal(vis=selfcal_vis, caltable='amp_6_longbaselines.cal', solint='30s', gaintype='G',
        uvrange='150~10000000m',
        calmode='a')

plotcal('amp_6_longbaselines.cal', xaxis='time', yaxis='amp', iteration='antenna',
        subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
        figfile='amp_6_longbaselines_vs_time.png')
plotcal('amp_6_longbaselines.cal', xaxis='antenna', yaxis='amp',
        figfile='amp_6_longbaselines_vs_antenna.png')

#applycal(vis=selfcal_vis, gaintable=["phase_6.cal", "amp_6_longbaselines.cal"],
#         interp="linear", applymode='calonly', calwt=False)







# dramatically overcleaned

# contimagename = 'Orion_SourceI_B6_continuum_r-2.clean0.4mJy.selfcal.ampphase6'
# os.system('rm -rf ' + contimagename + "*")
# tclean(vis=selfcal_vis,
#        imagename=contimagename,
#        field='Orion_BNKL_source_I',
#        specmode='mfs',
#        deconvolver='mtmfs',
#        nterms=2,
#        scales=[0,4,12],
#        smallscalebias=0.8,
#        imsize = imsize,
#        cell= cell,
#        weighting = 'briggs',
#        robust = -2,
#        niter = int(1e5),
#        threshold = '0.4mJy',
#        interactive = False,
#        outframe='LSRK',
#        veltype='radio',
#        savemodel='modelcolumn',
#        uvrange='50~36000m',
#       )
# makefits(contimagename)

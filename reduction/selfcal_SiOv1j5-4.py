
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



cell='0.004arcsec' # cell size for imaging.
imsize = [512,512] # size of image in pixels.

#contvis = ['band6_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
#if not os.path.exists(contvis[0]):
#    make_cont_ms()
#contvis = 'B6_calibrated_final_cont.ms'

restfreq_GHz = 215.59595

vis = 'band6.ms'
siovis = 'band6_siov1j5-4.ms'
if not os.path.exists(siovis):
    split(vis=vis, outputvis=siovis,
          spw='2:{0}~{1}GHz'.format(restfreq_GHz*(1-30/3e5),
                                    restfreq_GHz*(1+40/3e5),
                                   ),
          datacolumn='data'
         )

flagdata(vis=siovis, mode='manual', autocorr=True, uvrange='0~12m')


for imagename, depth, ii in [('Orion_SourceI_B6_SiO_v=1_J=5-4.noselfcal', '1Jy', 0),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase0', '300mJy', 1),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase1', '300mJy', 2),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase2', '200mJy', 3),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase3', '200mJy', 4),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase4', '100mJy', 5),
                             ('Orion_SourceI_B6_SiO_v=1_J=5-4.selfcal.phase5', '100mJy', 6),
                            ]:


    for suffix in ('', '.tt0', '.tt1', '.tt2'):
        for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
            todel = '{0}{1}{2}'.format(imagename, ext, suffix)
            if os.path.exists(todel):
                os.system('rm -rf {0}'.format(todel))


    tclean(vis=siovis,
           imagename=imagename,
           field='Orion_BNKL_source_I',
           specmode='cube',
           imsize = imsize,
           cell= cell,
           weighting = 'briggs',
           robust = -2,
           niter = int(1e5),
           threshold = '1.0Jy',
           interactive = False,
           outframe='LSRK',
           veltype='radio',
           savemodel='modelcolumn',
           uvrange='10~36000m',
           restfreq='{0}GHz'.format(restfreq_GHz),
           mask=['masercleancircle.crtf'],
          )

    makefits(imagename)


    rmtables(['sio_phase_{0}.cal'.format(ii)])
    gaincal(vis=siovis, caltable='sio_phase_{0}.cal'.format(ii), solint='int', gaintype='G',
            calmode='p')

    plotcal('sio_phase_{0}.cal'.format(ii), xaxis='time', yaxis='phase',
            iteration='antenna',
            subplot=331, timerange='2017/09/19/12:00:00~2017/09/20/12:00:00',
            figfile='sio_phase_{0}_vs_time.png'.format(ii))

    clearcal(vis=siovis)

    applycal(vis=siovis, gaintable=["phase_{0}.cal".format(ii)],
             interp="linear", applymode='calonly', calwt=False)

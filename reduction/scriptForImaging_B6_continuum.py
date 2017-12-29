
def make_cont_ms():

    mslist = ['uid___A001_X88e_X1d3_calibrated.ms',
              'uid___A002_Xb925ef_X4334_calibrated.ms']

    for ii,finalvis in enumerate(mslist):

        # Set spws to be used to form continuum
        contspws = '0,1,2,3'

        # If you have complex line emission and no dedicated continuum
        # windows, you will need to flag the line channels prior to averaging.
        flagmanager(vis=finalvis,mode='save',
                    versionname='before_cont_flags')

        initweights(vis=finalvis,wtmode='weight',dowtsp=True)

        # Flag the "line channels"
        flagchannels='0:80~110;155~195;235~260;325~355;420~440;480~495;585~610;695~720;770~790;865~885;1215~1230;1310~1495;1595~1615,1:0~170;315~450;585~600;670~715;845~890;960~985;1125~1875,2:35~90;305~345;470~490;530~605;660~780;865~1150;1260~1285;1450~1565;1685~1705;1785~1870,3:80~125;235~660;985~995;1045~1075;1235~1245;1585~1670;1765~1830'

        flagdata(vis=finalvis,mode='manual',
                  spw=flagchannels,flagbackup=False)


        # Average the channels within spws
        contvis='band6_continuum_ms{0}.ms'.format(ii)
        rmtables(contvis)
        os.system('rm -rf ' + contvis + '.flagversions')


        split2(vis=finalvis,
               spw=contspws,
               outputvis=contvis,
               width=[16,16,16,16], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
               datacolumn='data')


        # If you flagged any line channels, restore the previous flags
        flagmanager(vis=finalvis,mode='restore',
                    versionname='before_cont_flags')




def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)


if __name__ == "__main__":
    cell='0.004arcsec' # cell size for imaging.
    imsize = [7168,7168] # size of image in pixels.

    contvis = ['band6_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
    if not os.path.exists(contvis[0]):
        make_cont_ms()

    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B6_continuum_r{0}_dirty'.format(robust)

        for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
            rmtables(contimagename+ext)

        tclean(vis=contvis,
               imagename=contimagename,
               field='Orion_BNKL_source_I',
               specmode='mfs',
               deconvolver='mtmfs',
               nterms=2,
               scales=[0,4,12],
               imsize = imsize,
               cell= cell,
               weighting = 'briggs',
               robust = robust,
               niter = 0,
               threshold = '10Jy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
              )

        makefits(contimagename)


    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B6_continuum_r{0}'.format(robust)

        for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
            rmtables(contimagename+ext)

        tclean(vis=contvis,
               imagename=contimagename,
               field='Orion_BNKL_source_I',
               specmode='mfs',
               deconvolver='mtmfs',
               nterms=2,
               scales=[0,4,12],
               imsize = imsize,
               cell= cell,
               weighting = 'briggs',
               robust = robust,
               niter = int(1e5),
               threshold = '10mJy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
              )

        makefits(contimagename)

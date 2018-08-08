
def make_cont_ms():

    mslist = ['band7.ms', 'band7_lb.ms']

    for ii,finalvis in enumerate(mslist):

        # Set spws to be used to form continuum
        contspws = '25,27,29,31'

        # If you have complex line emission and no dedicated continuum
        # windows, you will need to flag the line channels prior to averaging.
        flagmanager(vis=finalvis,
                    mode='save',
                    versionname='before_cont_flags')

        initweights(vis=finalvis,wtmode='weight',dowtsp=True)

        # Flag the "line channels"
        flagchannels = '25:8~37;41~65;82~89;108~313;325~337;374~413;415~427;431~433;455~475;489~494;498~503;508~582;612~615;624~630;640~643;665~667;673~675;705~778;806~809;823~840;845~909;916~919;924~979;991~1002;1026~1041;1045~1061;1076~1080;1088~1153;1191~1198;1203~1218;1222~1228;1239~1343;1375~1382;1418~1448;1457~1462;1480~1501;1520~1534;1544~1548;1560~1607;1615~1617;1623~1630;1692~1919,27:0~12;32~65;114~189;279~280;298~331;340~343;371~372;380~381;393~648;662~679;708~712;724~746;749~754;770~792;819~868;881~965;978~997;1003~1010;1034~1037;1045~1048;1053~1059;1072~1077;1097~1100;1129~1135;1147~1162;1180~1209;1248~1342;1357~1363;1394~1425;1437~1473;1540~1545;1563~1566;1571~1612;1621~1632;1637~1653;1668~1683;1708~1711;1716~1721;1732~1749;1767~1782;1798~1800;1818~1832;1839~1841;1847~1855;1861~1919,29:0~27;51~55;83~90;124~131;229~280;310~314;327~383;431~448;474~483;488~525;535~586;596~621;634~637;641~645;658~660;684~686;693~695;697~736;744~832;890~896;912~920;937~987;995~1008;1082~1085;1098~1112;1126~1127;1141~1230;1253~1255;1258~1274;1307~1309;1387~1389;1402~1405;1426~1440;1493~1500;1521~1524;1539~1616;1634~1639;1664~1670;1687~1688;1707~1709;1726~1731;1757~1761;1784~1785;1788~1803;1818~1859;1894~1919,31:0~79;85~89;141~144;149~190;214~244;250~264;307~311;317~322;336~361;404~410;412~417;479~515;523~543;550~555;562~565;599~607;637~645;656~694;709~720;753~755;769~782;824~851;858~858;860~884;894~912;920~967;977~983;989~992;996~1010;1025~1307;1328~1342;1349~1351;1359~1366;1385~1402;1422~1503;1510~1520;1527~1556;1563~1590;1606~1631;1638~1655;1665~1668;1688~1711;1760~1770;1783~1786;1797~1815;1826~1829;1835~1919'

        flagdata(vis=finalvis,
                 mode='manual',
                 spw=flagchannels,
                 flagbackup=False, # no backup because we're manually backing up w/flagmanager
                )


        # Average the channels within spws
        contvis='band7_continuum_ms{0}.ms'.format(ii)
        rmtables(contvis)
        os.system('rm -rf ' + contvis + '.flagversions')


        split2(vis=finalvis,
               spw=contspws,
               outputvis=contvis,
               width=[16,16,16,16], # number of channels to average together.
               # The final channel width should be less than 125MHz in Bands 3,
               # 4, and 6 and 250MHz in Band 7.
               datacolumn='corrected')


        # If you flagged any line channels, restore the previous flags
        flagmanager(vis=finalvis,mode='restore',
                    versionname='before_cont_flags')

    concatvis = ['band7_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
    contvis = 'B7_calibrated_final_cont.ms'
    rmtables(contvis)
    os.system('rm -rf ' + contvis + '.flagversions')
    concat(vis=concatvis, concatvis=contvis)




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


if __name__ == "__main__":
    cell='0.004arcsec' # cell size for imaging.
    imsize = [6720,6720] # size of image in pixels.

    extensions = ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum','.alpha']

    #contvis = ['band7_continuum_ms{0}.ms'.format(ii) for ii in range(2)]
    #if not os.path.exists(contvis[0]):
    #    make_cont_ms()
    contvis = 'B7_calibrated_final_cont.ms'
    if not os.path.exists(contvis):
        make_cont_ms()
    else:
        print("Continuum measurement set {0} already exists.".format(contvis))

    redo = True
        

    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B7_continuum_r{0}_dirty'.format(robust)

        if os.path.exists(contimagename+".image.tt0.pbcor") and not redo:
            print("Skipping {0}".format(contimagename))
            continue
        elif redo:
            for suffix in ('', '.tt0', '.tt1', '.tt2'):
                for ext in extensions:
                    todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                    if os.path.exists(todel):
                        os.system('rm -rf {0}'.format(todel))

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
               uvrange='50~36000m',
              )

        makefits(contimagename)


    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B7_continuum_r{0}'.format(robust)

        if os.path.exists(contimagename+".image.tt0.pbcor") and not redo:
            print("Skipping {0}".format(contimagename))
            continue
        elif redo:
            for suffix in ('', '.tt0', '.tt1', '.tt2'):
                for ext in extensions:
                    todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                    if os.path.exists(todel):
                        os.system('rm -rf {0}'.format(todel))

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
               uvrange='50~36000m',
              )

        makefits(contimagename)



    # deeper clean

    dirtyimage = 'Orion_SourceI_B7_continuum_r-2_dirty.image.tt0'
    ia.open(dirtyimage)
    ia.calcmask(mask='"{0}" > 0.004'.format(dirtyimage), name='B7_clean_mask_5.0mJy')
    ia.close()
    makemask(mode='copy', inpimage=dirtyimage,
             inpmask=dirtyimage+":B7_clean_mask_5.0mJy", output='B7_clean_5.0mJy.mask',
             overwrite=True)
    exportfits('B7_clean_5.0mJy.mask', 'B7_clean_5.0mJy.mask.fits', dropdeg=True, overwrite=True)

    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B7_continuum_r{0}.mask5mJy.clean4mJy'.format(robust)

        if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
            if redo:
                for suffix in ('', '.tt0', '.tt1', '.tt2'):
                    for ext in extensions:
                        todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                        if os.path.exists(todel):
                            os.system('rm -rf {0}'.format(todel))


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
                   mask='B7_clean_5.0mJy.mask',
                   threshold = '4mJy',
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange='50~36000m',
                  )

            makefits(contimagename)


    # exclude all short baselines
    # 500m ~ 0.5"
    # 250m ~ 1"
    contimagename = 'Orion_SourceI_B7_continuum_r-2_longbaselines'

    if os.path.exists(contimagename+".image.tt0.pbcor") and redo:
        for suffix in ('', '.tt0', '.tt1', '.tt2'):
            for ext in extensions:
                todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                if os.path.exists(todel):
                    os.system('rm -rf {0}'.format(todel))

    if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
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
               robust = -2,
               niter = int(1e5),
               threshold = '5mJy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
               uvrange='250~36000m',
              )

    makefits(contimagename)


    contimagename = 'Orion_SourceI_B7_continuum_r0_250to2500m.mask5mjy'

    if os.path.exists(contimagename+".image.tt0.pbcor") and redo:
        for suffix in ('', '.tt0', '.tt1', '.tt2'):
            for ext in extensions:
                todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                if os.path.exists(todel):
                    os.system('rm -rf {0}'.format(todel))

    if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
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
               robust = 0,
               niter = int(1e5),
               threshold = '5mJy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
               uvrange='250~2500m',
               mask='B7_clean_5.0mJy.mask',
              )

    makefits(contimagename)


    contimagename = 'Orion_SourceI_B7_continuum_r0_250to5000m.mask5mjy'

    if os.path.exists(contimagename+".image.tt0.pbcor") and redo:
        for suffix in ('', '.tt0', '.tt1', '.tt2'):
            for ext in extensions:
                todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                if os.path.exists(todel):
                    os.system('rm -rf {0}'.format(todel))

    if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
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
               robust = 0,
               niter = int(1e5),
               threshold = '5mJy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
               uvrange='250~5000m',
               mask='B7_clean_5.0mJy.mask',
              )

    makefits(contimagename)


    for robust in (-2, 0, 2):
        contimagename = 'Orion_SourceI_B7_continuum_r{0}_1000to36000m.mask5mjy'.format(robust)

        if os.path.exists(contimagename+".image.tt0.pbcor") and redo:
            for suffix in ('', '.tt0', '.tt1', '.tt2'):
                for ext in extensions:
                    todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                    if os.path.exists(todel):
                        os.system('rm -rf {0}'.format(todel))

        if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
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
                   threshold = '5mJy',
                   interactive = False,
                   outframe='LSRK',
                   veltype='radio',
                   savemodel='none',
                   uvrange='1000~36000m',
                   mask='B7_clean_5.0mJy.mask',
                  )

        makefits(contimagename)


    # try auto-multithresh
    contimagename = 'Orion_SourceI_B7_continuum_r-2_automultithresh_1mJy'

    if os.path.exists(contimagename+".image.tt0.pbcor") and redo:
        for suffix in ('', '.tt0', '.tt1', '.tt2'):
            for ext in extensions:
                todel = '{0}{1}{2}'.format(contimagename, ext, suffix)
                if os.path.exists(todel):
                    os.system('rm -rf {0}'.format(todel))

    if redo or not os.path.exists(contimagename+".image.tt0.pbcor"):
        tclean(vis=contvis,
               imagename=contimagename,
               field='Orion_BNKL_source_I',
               specmode='mfs',
               deconvolver='mtmfs',
               nterms=2,
               scales=[0,4,12,36],
               imsize = imsize,
               cell= cell,
               weighting = 'briggs',
               robust = -2,
               niter = int(1e5),
               threshold = '1mJy',
               interactive = False,
               usemask='auto-multithresh',
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
               uvrange='50~36000m',
              )

    makefits(contimagename)

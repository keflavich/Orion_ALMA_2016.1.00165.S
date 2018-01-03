
def make_cont_ms():
    finalvis = 'member.uid___A001_X88e_X1d9_calibrated.ms'

    flagmanager(vis=finalvis, mode='save',
                versionname='before_cont_flags')

    initweights(vis=finalvis, wtmode='weight', dowtsp=True)

    # copied from scriptForImaging
    contspws = '0,1,2,3,4,5,6,7'
    flagchannels=('0:7~9;10~12;20~28;98~101;129~131;179~184;197~204;226~288;313~320;364~372;376~380;399~402;443~450;488~536;596~611;706~726;737~746;810~812;843~854;869~906;967~983;1013~1031;1094~1100;1113~1129;1138~1144;1154~1158;1207~1215;1259~1288;1305~1322;1334~1352;1391~1409;1429~1457;1587~1635;1657~1665;1682~1692;1709~1714;1722~1724;1740~1742;1805~1819;1825~1840;1871~1883,'
                  '1:12~59;92~99;101~124;147~155;225~249;269~273;299~303;341~362;378~382;388~394;398~403;472~476;481~495;521~543;558~569;603~630;651~658;708~716;735~748;764~779;844~863;893~899;916~919;925~941;1014~1024;1034~1041;1075~1103;1169~1192;1205~1211;1222~1226;1241~1274;1289~1294;1299~1308;1330~1347;1364~1373;1395~1403;1416~1422;1498~1505;1551~1567;1661~1684;1690~1695;1705~1711;1720~1726;1765~1771;1827~1919,'
                  '2:0~5;70~75;111~127;133~137;155~159;186~192;199~203;215~226;230~267;278~283;292~297;382~399;411~446;515~550;576~720;723~748;758~761;777~799;806~839;864~869;931~963;974~1000;1013~1148;1159~1283;1290~1324;1339~1387;1404~1412;1418~1481;1533~1560;1570~1610;1621~1630;1640~1652;1661~1665;1689~1699;1705~1713;1734~1810;1829~1838;1861~1894;1902~1908,'
                  '3:12~17;24~68;112~120;126~136;156~178;204~258;284~289;292~341;368~388;396~404;414~428;446~501;516~541;612~620;628~634;663~706;720~745;750~755;766~772;792~806;808~892;898~913;951~997;1008~1034;1080~1108;1121~1166;1275~1292;1303~1316;1326~1368;1377~1392;1405~1410;1422~1447;1535~1583;1620~1634;1665~1672;1678~1740;1768~1788;1803~1812,4:7~9;10~12;20~28;98~101;129~131;179~184;197~204;226~288;313~320;364~372;376~380;399~402;443~450;488~536;596~611;706~726;737~746;810~812;843~854;869~906;967~983;1013~1031;1094~1100;1113~1129;1138~1144;1154~1158;1207~1215;1259~1288;1305~1322;1334~1352;1391~1409;1429~1457;1587~1635;1657~1665;1682~1692;1709~1714;1722~1724;1740~1742;1805~1819;1825~1840;1871~1883, 5:12~59;92~99;101~124;147~155;225~249;269~273;299~303;341~362;378~382;388~394;398~403;472~476;481~495;521~543;558~569;603~630;651~658;708~716;735~748;764~779;844~863;893~899;916~919;925~941;1014~1024;1034~1041;1075~1103;1169~1192;1205~1211;1222~1226;1241~1274;1289~1294;1299~1308;1330~1347;1364~1373;1395~1403;1416~1422;1498~1505;1551~1567;1661~1684;1690~1695;1705~1711;1720~1726;1765~1771;1827~1919,6:0~5;70~75;111~127;133~137;155~159;186~192;199~203;215~226;230~267;278~283;292~297;382~399;411~446;515~550;576~720;723~748;758~761;777~799;806~839;864~869;931~963;974~1000;1013~1148;1159~1283;1290~1324;1339~1387;1404~1412;1418~1481;1533~1560;1570~1610;1621~1630;1640~1652;1661~1665;1689~1699;1705~1713;1734~1810;1829~1838;1861~1894;1902~1908,7:12~17;24~68;112~120;126~136;156~178;204~258;284~289;292~341;368~388;396~404;414~428;446~501;516~541;612~620;628~634;663~706;720~745;750~755;766~772;792~806;808~892;898~913;951~997;1008~1034;1080~1108;1121~1166;1275~1292;1303~1316;1326~1368;1377~1392;1405~1410;1422~1447;1535~1583;1620~1634;1665~1672;1678~1740;1768~1788;1803~1812'
                 )

    flagdata(vis=finalvis, mode='manual', spw=flagchannels, flagbackup=False)

    contvis='B3_calibrated_final_cont.ms'
    rmtables(contvis)
    os.system('rm -rf ' + contvis + '.flagversions')


    split2(vis=finalvis,
         spw=contspws,
         outputvis=contvis,
         width=[128,128,128,128,128,128,128,128], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
         datacolumn='data')


    flagmanager(vis=finalvis,mode='restore',
                versionname='before_cont_flags')


def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0', fitsimage=myimagebase+'.image.tt0.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.psf.tt0', fitsimage=myimagebase+'.psf.tt0.fits', dropdeg=True, overwrite=True) # export the PSF image


if __name__ == "__main__":
    # these comments copied from the QA2 delivery
    #206265/2.5e6 ~ 0.08
    #0.08/5 ~ 0.016 arcsec

    #6300/86 ~ 73
    #73/0.016 ~ 4608

    cell='0.016arcsec' # cell size for imaging.
    imsize = [4608,4608] # size of image in pixels.

    params = {2: {'imsize': [4608,4608], 'cell': '0.016arcsec'},
              0.5: {'imsize': [4608,4608], 'cell': '0.016arcsec'},
              -2: {'imsize': [4608*2,4608*2], 'cell': '0.008arcsec'},
             }



    contvis = 'B3_calibrated_final_cont.ms'
    if not os.path.exists(contvis):
        make_cont_ms()

    redo = False

    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B3_continuum_r{0}_dirty'.format(robust)

        imsize = params[robust]['imsize']
        cell = params[robust]['cell']

        if os.path.exists(contimagename+".image.tt0.pbcor") and not redo:
            continue
        elif redo:
            for suffix in ('', '.tt0', '.tt1', '.tt2'):
                for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
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
              )

        makefits(contimagename)


    for robust in (-2, 0.5, 2):
        contimagename = 'Orion_SourceI_B3_continuum_r{0}'.format(robust)

        imsize = params[robust]['imsize']
        cell = params[robust]['cell']

        if os.path.exists(contimagename+".image.tt0.pbcor") and not redo:
            continue

        elif redo:
            for suffix in ('', '.tt0', '.tt1', '.tt2'):
                for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
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
               threshold = '1mJy',
               interactive = False,
               outframe='LSRK',
               veltype='radio',
               savemodel='none',
              )

        makefits(contimagename)



    # deeper clean for robust -2 data
 
    dirtyimage = 'Orion_SourceI_B3_continuum_r-2_dirty.image.tt0'
    ia.open(dirtyimage)
    ia.calcmask(mask='"{0}" > 0.002'.format(dirtyimage), name='B3_clean_mask_2.0mJy')
    ia.close()
    makemask(mode='copy', inpimage=dirtyimage,
             inpmask=dirtyimage+":B3_clean_mask_2.0mJy", output='B3_clean_2.0mJy.mask',
             overwrite=True)
    exportfits('B3_clean_2.0mJy.mask', 'B3_clean_2.0mJy.mask.fits', dropdeg=True, overwrite=True)
 
    robust = -2
    contimagename = 'Orion_SourceI_B3_continuum_r{0}.mask2mJy.clean1mJy'.format(robust)
 
    imsize = params[robust]['imsize']
    cell = params[robust]['cell']
 
    for suffix in ('', '.tt0', '.tt1', '.tt2'):
        for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
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
           mask='B3_clean_2.0mJy.mask',
           threshold = '1mJy',
           interactive = False,
           outframe='LSRK',
           veltype='radio',
           savemodel='none',
         )

    makefits(contimagename)

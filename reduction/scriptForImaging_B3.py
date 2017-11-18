

########################################
# Check CASA version

import re
import casadef

if casadef.casa_version < '4.4.0' :
    sys.exit("Please use CASA version greater than or equal to 4.4.0 with this script")


##################################################
# Create an Averaged Continuum MS


finalvis='calibrated_final.ms' # This is your output ms from the data
                               # preparation script.

# Use plotms to identify line and continuum spectral windows.
plotms(vis=finalvis, xaxis='channel', yaxis='amplitude',
       ydatacolumn='data',
       avgtime='1e8', avgscan=True, avgchannel='1', 
       iteraxis='spw' )






# Set spws to be used to form continuum
contspws = '0,1,2,3,4,5,6,7'

# If you have complex line emission and no dedicated continuum
# windows, you will need to flag the line channels prior to averaging.
flagmanager(vis=finalvis,mode='save',
            versionname='before_cont_flags')

initweights(vis=finalvis,wtmode='weight',dowtsp=True)

# Flag the "line channels"
flagchannels='0:7~9;10~12;20~28;98~101;129~131;179~184;197~204;226~288;313~320;364~372;376~380;399~402;443~450;488~536;596~611;706~726;737~746;810~812;843~854;869~906;967~983;1013~1031;1094~1100;1113~1129;1138~1144;1154~1158;1207~1215;1259~1288;1305~1322;1334~1352;1391~1409;1429~1457;1587~1635;1657~1665;1682~1692;1709~1714;1722~1724;1740~1742;1805~1819;1825~1840;1871~1883, 1:12~59;92~99;101~124;147~155;225~249;269~273;299~303;341~362;378~382;388~394;398~403;472~476;481~495;521~543;558~569;603~630;651~658;708~716;735~748;764~779;844~863;893~899;916~919;925~941;1014~1024;1034~1041;1075~1103;1169~1192;1205~1211;1222~1226;1241~1274;1289~1294;1299~1308;1330~1347;1364~1373;1395~1403;1416~1422;1498~1505;1551~1567;1661~1684;1690~1695;1705~1711;1720~1726;1765~1771;1827~1919,2:0~5;70~75;111~127;133~137;155~159;186~192;199~203;215~226;230~267;278~283;292~297;382~399;411~446;515~550;576~720;723~748;758~761;777~799;806~839;864~869;931~963;974~1000;1013~1148;1159~1283;1290~1324;1339~1387;1404~1412;1418~1481;1533~1560;1570~1610;1621~1630;1640~1652;1661~1665;1689~1699;1705~1713;1734~1810;1829~1838;1861~1894;1902~1908,3:12~17;24~68;112~120;126~136;156~178;204~258;284~289;292~341;368~388;396~404;414~428;446~501;516~541;612~620;628~634;663~706;720~745;750~755;766~772;792~806;808~892;898~913;951~997;1008~1034;1080~1108;1121~1166;1275~1292;1303~1316;1326~1368;1377~1392;1405~1410;1422~1447;1535~1583;1620~1634;1665~1672;1678~1740;1768~1788;1803~1812,4:7~9;10~12;20~28;98~101;129~131;179~184;197~204;226~288;313~320;364~372;376~380;399~402;443~450;488~536;596~611;706~726;737~746;810~812;843~854;869~906;967~983;1013~1031;1094~1100;1113~1129;1138~1144;1154~1158;1207~1215;1259~1288;1305~1322;1334~1352;1391~1409;1429~1457;1587~1635;1657~1665;1682~1692;1709~1714;1722~1724;1740~1742;1805~1819;1825~1840;1871~1883, 5:12~59;92~99;101~124;147~155;225~249;269~273;299~303;341~362;378~382;388~394;398~403;472~476;481~495;521~543;558~569;603~630;651~658;708~716;735~748;764~779;844~863;893~899;916~919;925~941;1014~1024;1034~1041;1075~1103;1169~1192;1205~1211;1222~1226;1241~1274;1289~1294;1299~1308;1330~1347;1364~1373;1395~1403;1416~1422;1498~1505;1551~1567;1661~1684;1690~1695;1705~1711;1720~1726;1765~1771;1827~1919,6:0~5;70~75;111~127;133~137;155~159;186~192;199~203;215~226;230~267;278~283;292~297;382~399;411~446;515~550;576~720;723~748;758~761;777~799;806~839;864~869;931~963;974~1000;1013~1148;1159~1283;1290~1324;1339~1387;1404~1412;1418~1481;1533~1560;1570~1610;1621~1630;1640~1652;1661~1665;1689~1699;1705~1713;1734~1810;1829~1838;1861~1894;1902~1908,7:12~17;24~68;112~120;126~136;156~178;204~258;284~289;292~341;368~388;396~404;414~428;446~501;516~541;612~620;628~634;663~706;720~745;750~755;766~772;792~806;808~892;898~913;951~997;1008~1034;1080~1108;1121~1166;1275~1292;1303~1316;1326~1368;1377~1392;1405~1410;1422~1447;1535~1583;1620~1634;1665~1672;1678~1740;1768~1788;1803~1812' 

flagdata(vis=finalvis,mode='manual',
          spw=flagchannels,flagbackup=False)

# check that flags are as expected, NOTE must check reload on plotms
# gui if its still open.
plotms(vis=finalvis,yaxis='amp',xaxis='channel',
       avgchannel='1',avgtime='1e8',avgscan=True,iteraxis='spw') 

# Average the channels within spws
contvis='calibrated_final_cont.ms'
rmtables(contvis)
os.system('rm -rf ' + contvis + '.flagversions')


split2(vis=finalvis,
     spw=contspws,      
     outputvis=contvis,
     width=[128,128,128,128,128,128,128,128], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
     datacolumn='data')


# Check the weights. You will need to change antenna and field to
# appropriate values
plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DV24',field='3')

# If you flagged any line channels, restore the previous flags
flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

# Inspect continuum for any problems
plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')

# #############################################
# Image Parameters


# source parameters
# ------------------

field='3' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
imagermode='csclean' # uncomment if single field 
# imagermode='mosaic' # uncomment if mosaic or if combining one 7m and one 12m pointing.
# phasecenter=3 # uncomment and set to field number for phase
                # center. Note lack of ''.  Use the weblog to
                # determine which pointing to use. Remember that the
                # field ids for each pointing will be re-numbered
                # after your initial split. You can also specify the
                # phase center using coordinates, e.g.,
                # phasecenter='J2000 19h30m00 -40d00m00'

# image parameters.
# ----------------





#206265/2.5e6 ~ 0.08
#0.08/5 ~ 0.016 arcsec

#6300/86 ~ 73
#73/0.016 ~ 4608

cell='0.016arcsec' # cell size for imaging.
imsize = [4608,4608] # size of image in pixels.

# velocity parameters
# -------------------

outframe='LSRK' # velocity reference frame. See science goals.
veltype='radio' # velocity type. 


# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'

#############################################
# Imaging the Continuuum

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'         
contimagename = 'calibrated_final_cont'

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename+ext)



#Clean Cycles: 4
#Beam Area: 0.10 by 0.07 arcsec
#RMS:  20 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)



# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an optional step.
#contmaskname = 'cont.mask'
###rmtables(contmaskname) # if you want to delete the old mask
#os.system('cp -ir ' + contimagename + '.mask ' + contmaskname)

##############################################
# Self-calibration on the continuum [OPTIONAL]



contvis = 'calibrated_final_cont.ms'         
contimagename = 'calibrated_final_cont'

refant = 'DV25' # reference antenna.


spwmap = [0,0,0,0,0,0,0,0] # mapping self-calibration solutions to individual spectral windows. Generally an array of n zeroes, where n is the number of spectral windows in the data sets.

# save initial flags in case you don't like the final
# self-calibration. The task applycal will flag data that doesn't have
# solutions.
flagmanager(vis=contvis,mode='save',versionname='before_selfcal',merge='replace')

# Get rid of any models that might be hanging around in the image header
delmod(vis=contvis,field=field,otf=True)

# If you are re-doing your self-cal, uncomment the next line to reset
# your corrected data column back to its original state.
#clearcal(vis=contvis)

# shallow clean on the continuum

for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename + '_p0'+ ext)
    
#Clean Cycles: 1
#Beam Area: 0.10 by 0.07 arcsec
#RMS: 25 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      imagename=contimagename + '_p0',
      field=field,
#      phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)


# per scan solution
rmtables('pcal1')
gaincal(vis=contvis,
        caltable='pcal1',
        field=field,
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='inf',
        minsnr=3.0,
        minblperant=6)

# Check the solution
plotcal(caltable='pcal1',
        xaxis='time',
        yaxis='phase',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,-180,180])

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         field=field,
         spwmap=spwmap, 
         gaintable=['pcal1'],
         gainfield='',
         calwt=False, 
         flagbackup=False,
         interp='linearperobs')

# clean deeper
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename + '_p1'+ ext)

#Clean Cycles: 3
#Beam Area: 0.10 by 0.07 arcsec
#RMS: 21 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      field=field,
#      phasecenter=phasecenter, # uncomment if mosaic.      
      imagename=contimagename + '_p1',
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)

# Note number of iterations performed.

# shorter solution
rmtables('pcal2')
gaincal(vis=contvis,
        field=field,
        caltable='pcal2',
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='30.25s', # solint=30.25s gets you five 12m integrations, while solint=50.5s gets you five 7m integration
        minsnr=3.0,
        minblperant=6)

# Check the solution
plotcal(caltable='pcal2',
        xaxis='time',
        yaxis='phase',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,-180,180])

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         spwmap=spwmap, 
         field=field,
         gaintable=['pcal2'],
         gainfield='',
         calwt=False, 
         flagbackup=False,
         interp='linearperobs')

# clean deeper
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename + '_p2'+ ext)

#Clean Cycles: 4 
#Beam Area: 0.10 by 0.07 arcsec
#RMS: 21 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      imagename=contimagename + '_p2',
      field=field,
#      phasecenter=phasecenter, # uncomment if mosaic.            
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)


# shorter solution
rmtables('pcal3')
gaincal(vis=contvis,
        field=field,
        caltable='pcal3',
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='int',
        minsnr=3.0,
        minblperant=6)

# Check the solution
plotcal(caltable='pcal3',
        xaxis='time',
        yaxis='phase',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,-180,180])

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         spwmap=spwmap,
         field=field,
         gaintable=['pcal3'],
         gainfield='',
         calwt=False, 
         flagbackup=False,
         interp='linearperobs')

# do the amplitude self-calibration.
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename + '_p3'+ ext)

#Clean Cycles: 6
#Beam Area: 0.10 by 0.07 arcsec
#RMS: 20 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      imagename=contimagename + '_p3',
      field=field,
#      phasecenter=phasecenter, # uncomment if mosaic.            
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)


rmtables('apcal')
gaincal(vis=contvis,
        field=field,
        caltable='apcal',
        gaintype='T',
        refant=refant,
        calmode='ap',
        combine='spw',
        solint='inf',
        minsnr=3.0,
        minblperant=6,
#        uvrange='>50m', # may need to use to exclude extended emission
        gaintable='pcal3',
        spwmap=spwmap,
        solnorm=True)

plotcal(caltable='apcal',
        xaxis='time',
        yaxis='amp',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,0.2,1.8])

applycal(vis=contvis,
         spwmap=[spwmap,spwmap], # select which spws to apply the solutions for each table
         field=field,
         gaintable=['pcal3','apcal'],
         gainfield='',
         calwt=False,
         flagbackup=False,
         interp=['linearperobs','linearperobs'])

# Make amplitude and phase self-calibrated image.
for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
    rmtables(contimagename + '_ap'+ ext)

#Clean Cycles: 4
#Beam Area: 0.10 by 0.07 arcsec
#RMS: 18 uJy/Beam for 3.36 GHz Bandwidth
clean(vis=contvis,
      imagename=contimagename + '_ap',
      field=field, 
#      phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)


# Save results of self-cal in a new ms
split(vis=contvis,
      outputvis=contvis+'.selfcal',
      datacolumn='corrected')

# reset the corrected data column in the  ms to the original calibration.

clearcal(vis=contvis)



#########################################
## Continuum Subtraction for Line Imaging
#
##>>> If you have observations that include both line and strong (>3 sigma
##>>> per final line image channel) continuum emission, you need to
##>>> subtract the continuum from the line data. You should not continuum
##>>> subtract if the line of interest is in absorption.
#
##>>> NOTE THAT WE'RE BACK TO EXCLUDECHANS=FALSE, SO FITSPW INDICATES THE LINE-FREE CHANNELS.
#
##>>> You can use au.invertChannelRanges(flagchannels,vis=finalvis) to
##>>> get the fitspw below. You will need to insert any continuum spws
##>>> that weren't included in flagchannels. For example, if your continuum
##>>> spws are '0,1,2' and flagchannels='1:260~500', au.invertChannelRanges will return
##>>> '1:0~259,1:501~3839'. The fitspw parameter should be '0,1:0~259,1:501~3839,2'
##>>> Make sure to cut and paste the output in fitspw below since PIs don't have
##>>> analysisUtilities by default.
#
#fitspw = '0,1,2:0~1200;1500~3839,3:0~1200;1500~3839' # *line-free* channels for fitting continuum
#linespw = '2,3' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once.
#
#finalvis='calibrated_final.ms'
#
#uvcontsub(vis=finalvis,
#          spw=linespw, # spw to do continuum subtraction on
#          fitspw=fitspw, # regions without lines.
#          excludechans=False, # fit the regions in fitspw
#          combine='spw', 
#          solint='int',
#          fitorder=1,
#          want_cont=False) # This value should not be changed.
#
##>>> Note that the continuum subtraction is done for each field in 
##>>> turn. However, if the fields have different line-free channels, you
##>>> will need to do the continuum subtraction separately for each field.
#
## NOTE: Imaging the continuum produced by uvcontsub with
## want_cont=True will lead to extremely poor continuum images because
## of bandwidth smearing effects. For imaging the continuum, you should
## always create a line-free continuum data set using the process
## outlined above.
#
##########################################################
## Apply continuum self-calibration to line data [OPTIONAL]
#
## uncomment one  of the following
## linevis = finalvis+'.contsub' # if continuum subtracted
## linevis = finalvis  #  if not continuum subtracted
#
## save original flags in case you don't like the self-cal
#flagmanager(vis=linevis,mode='save',versionname='before_selfcal',merge='replace')
#
#spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
#applycal(vis=linevis,
#         spwmap=[spwmap_line, spwmap_line], # entering the appropriate spwmap_line value for each spw in the input dataset
#         field=field,
#         gaintable=['pcal3','apcal'],
#         gainfield='',
#         calwt=False,
#         flagbackup=False,
#         interp=['linearperobs','linearperobs'])
#
## Save results of self-cal in a new ms and reset the image name.
#split(vis=linevis,
#      outputvis=linevis+'.selfcal',
#      datacolumn='corrected')
#
## reset the corrected data column in the  ms to the original calibration
##>>> This can also be used to return your ms to it's original
##>>> pre-self-cal state if you are unhappy with your self-calibration.
#clearcal(linevis)
#
##>>> The applycal task will automatically flag data without good
##>>> gaincal solutions. If you are unhappy with your self-cal and wish to
##>>> return the flags to their original state, run the following command
##>>> flagmanager(vis=linevis, mode='restore',versionname='before_selfcal')
#
#linevis=linevis+'.selfcal'
#
###############################################
## Image line emission [REPEAT AS NECESSARY]
#
##>>> If you did an mstransform/cvel, use the same velocity parameters in
##>>> the clean that you did for the regridding. If you did not do an
##>>> mstransform and have multiple executions of a scheduling block,
##>>> select the spws with the same rest frequency using the spw parameter
##>>> (currently commented out below). DO NOT INCLUDE SPWS WITH DIFFERENT
##>>> REST FREQUENCIES IN THE SAME RUN OF CLEAN: THEY WILL SLOW DOWN
##>>> IMAGING CONSIDERABLY.
#
#finalvis = 'calibrated_final.ms'
## linevis = finalvis # uncomment if you neither continuum subtracted nor self-calibrated your data.
## linevis = finalvis + '.contsub' # uncomment if continuum subtracted
## linevis = finalvis + '.contsub.selfcal' # uncommment if both continuum subtracted and self-calibrated
## linevis = finalvis + '.selfcal' # uncomment if just self-calibrated (no continuum subtraction)
#
#sourcename ='n253' # name of source
#linename = 'CO10' # name of transition (see science goals in OT for name) 
#lineimagename = sourcename+'_'+linename # name of line image
#
#restfreq='115.27120GHz' # Typically the rest frequency of the line of
#                        # interest. If the source has a significant
#                        # redshift (z>0.2), use the observed sky
#                        # frequency (nu_rest/(1+z)) instead of the
#                        # rest frequency of the
#                        # line.
#
## spw='1' # uncomment and replace with appropriate spw if necessary.
#
#start='-100km/s' # start velocity. See science goals for appropriate value.
#width='2km/s' # velocity width. See science goals.
#nchan = 100  # number of channels. See science goals for appropriate value.
#
##>>> To specify a spws from multiple executions that had not been regridded using cvel, use
##>>>       import numpy as np
##>>>       spw = str.join(',',map(str,np.arange(0,n,nspw)))
##>>>
##>>> where n is the total number of windows x executions and nspw is the
##>>> number of spectral windows per execution. Note that the spectral
##>>> windows need to have the same order in all data sets for this code
##>>> to work. Add a constant offset (i.e., +1,+2,+3) to the array
##>>> generated by np.arange to get the other sets of windows.
#
## If necessary, run the following commands to get rid of older clean
## data.
#
##clearcal(vis=linevis)
##delmod(vis=linevis)
#
#for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#    rmtables(lineimagename + ext)
#
#clean(vis=linevis,
#      imagename=lineimagename, 
#      field=field,
##      spw=spw,
##      phasecenter=phasecenter, # uncomment if mosaic.      
#      mode='velocity',
#      start=start,
#      width=width,
#      nchan=nchan, 
#      outframe=outframe, 
#      veltype=veltype, 
#      restfreq=restfreq, 
#      niter=niter,  
#      threshold=threshold, 
#      interactive=True,
#      cell=cell,
#      imsize=imsize, 
#      weighting=weighting, 
#      robust=robust,
#      imagermode=imagermode)
#
##>>> If interactively cleaning (interactive=True), then note number of
##>>> iterations at which you stop for the PI. This number will help the
##>>> PI replicate the delivered images. Do not clean empty
##>>> images. Just click the red X to stop the interactive and note the
##>>> RMS.
#
## If you'd like to redo your clean, but don't want to make a new mask
## use the following commands to save your original mask. This is an
## optional step.
## linemaskname = 'line.mask'
### rmtables(linemaskname) # uncomment if you want to overwrite the mask.
## os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)
#
##############################################
# Apply a primary beam correction

import glob

myimages = glob.glob("*.image")

rmtables('*.pbcor')
for image in myimages:
    pbimage = image.rsplit('.',1)[0]+'.flux'
    outfile = image.rsplit('.',1)[0]+'.pbcor'
    impbcor(imagename=image, pbimage=pbimage, outfile = outfile)

##############################################
# Export the images

import glob

myimages = glob.glob("*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*.flux")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs


os.system("rm -rf *.png")
mycontimages = glob.glob("calibrated*.image")
for cimage in mycontimages:
    mymax=imstat(cimage)['max'][0]
    mymin=-0.1*mymax
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage,'range':[mymin,mymax]},out=outimage)


# this will have to be run for each sourcename
#sourcename='' # insert source here, if it isn't already set
#mylineimages = glob.glob(sourcename+"*.image")
#for limage in mylineimages:
#    rms=imstat(limage,chans='1')['rms'][0]
#    mom8=limage+'.mom8'
#    os.system("rm -rf "+mom8)
#    immoments(limage,moments=[8],outfile=mom8)
#    mymax=imstat(mom8)['max'][0]
#    mymin=-0.1*mymax
#    os.system("rm "+mom8+".png")
#    imview(raster={'file':mom8,'range':[mymin,mymax]},out=mom8+'.png')


##############################################
# Analysis

# For examples of how to get started analyzing your data, see
#     https://casaguides.nrao.edu/index.php/TWHydraBand7_Imaging_4.3
#     

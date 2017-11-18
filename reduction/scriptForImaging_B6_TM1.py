


########################################
# Check CASA version

import re

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
     width=[16,16,16,16], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
     datacolumn='data')


# Check the weights. You will need to change antenna and field to
# appropriate values
plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DV24',field='4')

# If you flagged any line channels, restore the previous flags
flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

# Inspect continuum for any problems
plotms(vis=contvis,xaxis='uvdist',yaxis='amp',coloraxis='spw')

# #############################################
# Image Parameters


# source parameters
# ------------------

field='4' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
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




cell='0.004arcsec' # cell size for imaging.
imsize = [7168,7168] # size of image in pixels.

# velocity parameters
# -------------------

outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. 


# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=-2.0
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



clean(vis=contvis,
      imagename=contimagename,
      field=field,
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
##rmtables(contmaskname) # if you want to delete the old mask
#os.system('cp -ir ' + contimagename + '.mask ' + contmaskname)

##############################################
# Self-calibration on the continuum [OPTIONAL]



# contvis = 'calibrated_final_cont.ms'         
# contimagename = 'calibrated_final_cont'

# refant = 'DV24' # reference antenna.

# #>>> Choose a reference antenna that's in the array. The tasks plotants
# #>>> and listobs/vishead can tell you what antennas are in the array. For
# #>>> data sets with multiple executions, you will want to choose an antenna
# #>>> that's present in all the executions. The task au.commonAntennas()
# #>>> can help with this.

# spwmap = [0,0,0,0] # mapping self-calibration solutions to individual spectral windows. Generally an array of n zeroes, where n is the number of spectral windows in the data sets.

# # save initial flags in case you don't like the final
# # self-calibration. The task applycal will flag data that doesn't have
# # solutions.
# flagmanager(vis=contvis,mode='save',versionname='before_selfcal',merge='replace')

# # Get rid of any models that might be hanging around in the image header
# delmod(vis=contvis,field=field,otf=True)

# # If you are re-doing your self-cal, uncomment the next line to reset
# # your corrected data column back to its original state.
# #clearcal(vis=contvis)

# # shallow clean on the continuum

# for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#     rmtables(contimagename + '_p0'+ ext)
    
# clean(vis=contvis,
#       imagename=contimagename + '_p0',
#       field=field,
#       mode='mfs',
#       psfmode='clark',
#       imsize = imsize, 
#       cell= cell, 
#       weighting = weighting, 
#       robust=robust,
#       niter=niter, 
#       threshold=threshold, 
#       interactive=True,
#       usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
#       imagermode=imagermode)



# #>>> Note number of iterations performed.

# # per scan solution
# rmtables('pcal1')
# gaincal(vis=contvis,
#         caltable='pcal1',
#         field=field,
#         gaintype='T',
#         refant=refant, 
#         calmode='p',
#         combine='spw', 
#         solint='inf',
#         minsnr=3.0,
#         minblperant=6)

# # Check the solution
# plotcal(caltable='pcal1',
#         xaxis='time',
#         yaxis='phase',
#         timerange='',
#         iteration='antenna',
#         subplot=421,
#         plotrange=[0,0,-180,180])

# # apply the calibration to the data for next round of imaging
# applycal(vis=contvis,
#          field=field,
#          spwmap=spwmap, 
#          gaintable=['pcal1'],
#          gainfield='',
#          calwt=False, 
#          flagbackup=False,
#          interp='linearperobs')

# # clean deeper
# for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#     rmtables(contimagename + '_p1'+ ext)

# clean(vis=contvis,
#       field=field,
#       imagename=contimagename + '_p1',
#       mode='mfs',
#       psfmode='clark',
#       imsize = imsize, 
#       cell= cell, 
#       weighting = weighting, 
#       robust=robust,
#       niter=niter, 
#       threshold=threshold, 
#       interactive=True,
#       usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
#       imagermode=imagermode)

# Note number of iterations performed.

# 200 iterations
# 0.05/0.00015 = 333

# shorter solution
# rmtables('pcal2')
# gaincal(vis=contvis,
#         field=field,
#         caltable='pcal2',
#         gaintype='T',
#         refant=refant, 
#         calmode='p',
#         combine='spw', 
#         solint='30.25s', # solint=30.25s gets you five 12m integrations, while solint=50.5s gets you five 7m integration
#         minsnr=3.0,
#         minblperant=6)

# # Check the solution
# plotcal(caltable='pcal2',
#         xaxis='time',
#         yaxis='phase',
#         timerange='',
#         iteration='antenna',
#         subplot=421,
#         plotrange=[0,0,-180,180])

# # apply the calibration to the data for next round of imaging
# applycal(vis=contvis,
#          spwmap=spwmap, 
#          field=field,
#          gaintable=['pcal2'],
#          gainfield='',
#          calwt=False, 
#          flagbackup=False,
#          interp='linearperobs')

# # clean deeper
# for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#     rmtables(contimagename + '_p2'+ ext)

# clean(vis=contvis,
#       imagename=contimagename + '_p2',
#       field=field,
#       mode='mfs',
#       psfmode='clark',
#       imsize = imsize, 
#       cell= cell, 
#       weighting = weighting, 
#       robust=robust,
#       niter=niter, 
#       threshold=threshold, 
#       interactive=True,
#       usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
#       imagermode=imagermode)

# 0.0325/0.0001 = 325

# shorter solution
# rmtables('pcal3')
# gaincal(vis=contvis,
#         field=field,
#         caltable='pcal3',
#         gaintype='T',
#         refant=refant, 
#         calmode='p',
#         combine='spw', 
#         solint='int',
#         minsnr=3.0,
#         minblperant=6)

# # Check the solution
# plotcal(caltable='pcal3',
#         xaxis='time',
#         yaxis='phase',
#         timerange='',
#         iteration='antenna',
#         subplot=421,
#         plotrange=[0,0,-180,180])

# # apply the calibration to the data for next round of imaging
# applycal(vis=contvis,
#          spwmap=spwmap,
#          field=field,
#          gaintable=['pcal3'],
#          gainfield='',
#          calwt=False, 
#          flagbackup=False,
#          interp='linearperobs')

# # do the amplitude self-calibration.
# for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#     rmtables(contimagename + '_p3'+ ext)

# clean(vis=contvis,
#       imagename=contimagename + '_p3',
#       field=field,
# #      phasecenter=phasecenter, # uncomment if mosaic.            
#       mode='mfs',
#       psfmode='clark',
#       imsize = imsize, 
#       cell= cell, 
#       weighting = weighting, 
#       robust=robust,
#       niter=niter, 
#       threshold=threshold, 
#       interactive=True,
#       usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
#       imagermode=imagermode)

# #>>> Note number of iterations performed.

# rmtables('apcal')
# gaincal(vis=contvis,
#         field=field,
#         caltable='apcal',
#         gaintype='T',
#         refant=refant,
#         calmode='ap',
#         combine='spw',
#         solint='inf',
#         minsnr=3.0,
#         minblperant=6,
# #        uvrange='>50m', # may need to use to exclude extended emission
#         gaintable='pcal3',
#         spwmap=spwmap,
#         solnorm=True)

# plotcal(caltable='apcal',
#         xaxis='time',
#         yaxis='amp',
#         timerange='',
#         iteration='antenna',
#         subplot=421,
#         plotrange=[0,0,0.2,1.8])

# applycal(vis=contvis,
#          spwmap=[spwmap,spwmap], # select which spws to apply the solutions for each table
#          field=field,
#          gaintable=['pcal3','apcal'],
#          gainfield='',
#          calwt=False,
#          flagbackup=False,
#          interp=['linearperobs','linearperobs'])

# # Make amplitude and phase self-calibrated image.
# for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage','.pb','.wtsum']:
#     rmtables(contimagename + '_ap'+ ext)

# clean(vis=contvis,
#       imagename=contimagename + '_ap',
#       field=field, 
# #      phasecenter=phasecenter, # uncomment if mosaic.      
#       mode='mfs',
#       psfmode='clark',
#       imsize = imsize, 
#       cell= cell, 
#       weighting = weighting, 
#       robust=robust,
#       niter=niter, 
#       threshold=threshold, 
#       interactive=True,
#       usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
#       imagermode=imagermode)


## Save results of self-cal in a new ms
#split(vis=contvis,
#      outputvis=contvis+'.selfcal',
#      datacolumn='corrected')

# reset the corrected data column in the  ms to the original calibration.

#clearcal(vis=contvis)





##############################################

# Apply a primary beam correction

import glob

myimages = ['calibrated_final_cont.image']

rmtables('*.pbcor')
for image in myimages:
    pbimage = image.rsplit('.',1)[0]+'.flux'
    outfile = image.rsplit('.',1)[0]+'.pbcor'
    impbcor(imagename=image, pbimage=pbimage, outfile = outfile)

##############################################
# Export the images

import glob

myimages = ["calibrated_final_cont"]
for image in myimages:
    exportfits(imagename=image+'.pbcor', fitsimage=image+'.pbcor.fits',overwrite=True)

for image in myimages:
    exportfits(imagename=image+'.flux', fitsimage=image+'.flux.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs


os.system("rm -rf *.png")
mycontimages = glob.glob("calibrated*.image")
for cimage in mycontimages:
    max=imstat(cimage)['max'][0]
    min=-0.1*max
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage,'range':[min,max]},out=outimage)


##############################################
# Analysis

# For examples of how to get started analyzing your data, see
#     https://casaguides.nrao.edu/index.php/TWHydraBand7_Imaging_4.3

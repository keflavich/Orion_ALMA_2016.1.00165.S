
########################################
# Getting a list of ms files to image

import glob


##################################################
# Flag Bad Data [OPTIONAL]


# Save original flags
flagmanager(vis='uid___A002_Xb925ef_X4334_target.ms',
                mode='save',
                versionname='original_flags')

# Inspect the science data
#fieldlist = ['3'] # list of science data fields to inspect
#spwlist = ['1'] # list of science spws to inspect

# loop through science data fields and spws to inspect.

#for vis in vislist:
#    for field in fieldlist:
#        for spw in spwlist:
#            plotms(vis=vis,xaxis='uvwave',yaxis='amp',avgtime='3e8',
#                   field=field,spw=spw) 
#            raw_input("push enter to continue")
#
#            plotms(vis=vis,xaxis='chan',yaxis='amp',avgtime='3e8',
#                   field=field,spw=spw) 
#            raw_input("push enter to continue")

# Flag the offending data. See flagdata help for more info.
#flagdata(vis='',mode='manual',action='apply',flagbackup=False)

# If you need to restore original flags, use the following command.
#flagmanager(vis='',mode='restore',versionname='original_flags')

########################################
# Flux Equalization [OPTIONAL]


# generating the script -- REMOVE BEFORE SENDING TO PI
#es.generateReducScript(['uid_FIRST-EB.ms.split.cal','uid_SECOND-EB.ms.split.cal',(etc)], step='fluxcal')


###############################################################
# Combining Measurement Sets from Multiple Executions 


# If you have multiple executions, you will want to combine the
# scheduling blocks into a single ms using concat for ease of imaging
# and self-calibration. Each execution of the scheduling block will
# generate multiple spectral windows with different sky frequencies,
# but the same rest frequency, due to the motion of the Earth. Thus,
# the resulting concatentated file will contain n spws, where n is
# (#original science spws) x (number executions).  In other words, the
# multiple spws associated with a single rest frequency will not be
# regridded to a single spectral window in the ms.

#concatvis='calibrated.ms'

#rmtables(concatvis)
#os.system('rm -rf ' + concatvis + '.flagversions')
#concat(vis=vislist,
#       concatvis=concatvis)

###################################
# Splitting off science target data

# concatvis = vislist[0]

# concatvis='calibrated.ms'


#sourcevis='calibrated_source.ms'
#rmtables(sourcevis)
#os.system('rm -rf ' + sourcevis + '.flagversions')
#split(vis=concatvis,
#      intent='*TARGET*', # split off the target sources
#      outputvis=sourcevis,
#      datacolumn='data')

###############################################################
# Regridding spectral windows [OPTIONAL]

#sourcevis='calibrated_source.ms'
#regridvis='calibrated_source_regrid.ms'
#veltype = 'radio' # Keep set to radio. See notes in imaging section.
#width = '0.23km/s' # see science goals in the OT
#nchan = -1 # leave this as the default
#mode='velocity' # see science goals in the OT
#start='' # leave this as the default
#outframe = 'bary' # velocity reference frame. see science goals in the OT.
#restfreq='115.27120GHz' # rest frequency of primary line of interest. 
#field = '4' # select science fields.
#spw = '0,5,10' # spws associated with a single rest frequency. Do not attempt to combine spectral windows associated with different rest frequencies. This will take a long time to regrid and most likely isn't what you want.

#rmtables(regridvis)
#os.system('rm -rf ' + regridvis + '.flagversions')
    
#cvel(vis=sourcevis,
#     field=field,
#     outputvis=regridvis,
#     spw=spw,
#     mode=mode,
#     nchan=nchan,
#     width=width,
#     start=start,
#     restfreq=restfreq,
#     outframe=outframe,
#     veltype=veltype)


############################################
# Rename and backup data set

os.system('mv -i uid___A002_Xb925ef_X4334_target.ms' + ' ' + 'calibrated_final.ms')

# os.system('mv -i ' + regridvis + ' ' + 'calibrated_final.ms') 

# At this point you should create a backup of your final data set in
# case the ms you are working with gets corrupted by clean. 

os.system('cp -ir calibrated_final.ms calibrated_final.ms.backup')


############################################
# Output a listobs file

listobs(vis='calibrated_final.ms',listfile='calibrated_final.ms.listobs.txt') 

"""
This is a CASA imaging script, so it should be run from within casa with
%run -i or similar.  Therefore, CASA tasks are not imported.
"""
import numpy as np
import os
import glob
import datetime
from astropy import constants
from astropy import units as u

from line_to_image_list import mses, line_to_image_list, imaging_parameters, ms_basepath

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

    if cleanup:
        for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                       'mask', 'image', 'workdirectory'):
            os.system('rm -rf "{0}.{1}"'.format(myimagebase, suffix))


for line_info in line_to_image_list:

    v0 = u.Quantity(line_info['velocity_range'][0], u.km/u.s)
    v1 = u.Quantity(line_info['velocity_range'][1], u.km/u.s)
    band = 'b{0}'.format(line_info['band'])
    linename = line_info['name']
    frequency = line_info['frequency'].to(u.Hz).value

    # only use first MS for metadata, but use both later
    msfile = os.path.join(ms_basepath,
                          mses[band][0])

    msmd.open(msfile)
    spws = msmd.spwsforfield('Orion_BNKL_source_I')
    spws = [x for x in spws
            if (msmd.chanfreqs(x).min() < frequency) and
            (msmd.chanfreqs(x).max() > frequency)]
    dnu = np.diff(msmd.chanfreqs(spws[0])).mean()
    dv = (dnu / frequency * constants.c).to(u.km/u.s)
    nchan = int(np.abs(((v1-v0)/dv).decompose()))
    msmd.close()

    spw = ",".join([str(spw) for spw in spws])

    start_vel = v0
    print("Start_vel = {0}".format(start_vel))

    for suffix, niter in (('clarkclean1000', 1000), ):


        imagename = 'OrionFullField.{3}.spw{0}.{2}.{1}'.format(spw, suffix, linename, band)
        if not os.path.exists("{0}.image.pbcor.fits".format(imagename)):
            print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
            tclean(vis=mses[band],
                   imagename=imagename,
                   field='Orion_BNKL_source_I',
                   spw=spw,
                   gridder='standard',
                   specmode='cube',
                   start='{0}km/s'.format(start_vel.to(u.km/u.s).value),
                   nchan=nchan,
                   restfreq=frequency,
                   veltype='radio',
                   outframe='LSRK',
                   interactive=False,
                   niter=imaging_parameters[band]['niter'],
                   imsize=imaging_parameters[band]['imsize'],
                   cell=imaging_parameters[band]['cell'],
                   weighting='briggs',
                   robust=0.5,
                   phasecenter='',
                   threshold=imaging_parameters[band]['threshold'],
                   savemodel='none',
                  )
            makefits(imagename)

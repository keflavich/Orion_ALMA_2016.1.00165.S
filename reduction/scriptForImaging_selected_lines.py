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
    if line_info['name'] != 'SiS_12-11':
        continue

    v0 = u.Quantity(line_info['velocity_range'][0], u.km/u.s)
    v1 = u.Quantity(line_info['velocity_range'][1], u.km/u.s)
    band = 'b{0}'.format(line_info['band'])
    linename = line_info['name']
    frequency = line_info['frequency'].to(u.Hz).value

    # only use first MS for metadata, but use both later

    dvs = []
    spws = {}

    for ms_ in mses[band]:
        msfile = os.path.join(ms_basepath,
                              ms_)

        msmd.open(msfile)
        spws_ = msmd.spwsforfield('Orion_BNKL_source_I')
        spws_ = [x for x in spws_
                if (msmd.chanfreqs(x).min() < frequency) and
                (msmd.chanfreqs(x).max() > frequency)]
        if len(spws_) == 0:
            print("SKIPPING {0} BECAUSE NO SPWS FOUND!!".format(line_info))
            continue
        for spw in tuple(spws_):
            dnu = np.diff(msmd.chanfreqs(spw)).mean()
            if np.abs(dnu) > 2e6:
                spws_.remove(spw)
                continue
            dv = (dnu / frequency * constants.c).to(u.km/u.s)
            dvs.append(dv)
        dv = min(dvs)
        spws[ms_] = spws_

    width = '{0}km/s'.format(dv)

    nchan = int(np.abs(((v1-v0)/dv).decompose()))
    msmd.close()

    print(spws)

    spw = [",".join([str(spw) for spw in spwl]) for _,spwl in spws.items()]

    start_vel = v0
    print("Start_vel = {0} for line {1} restfreq={2}".format(start_vel, linename, frequency))
    print("nchan = {0}".format(nchan))
    print("spw = {0}".format(spw))
    spwname = (str(spw).replace("[","")
               .replace(",","_")
               .replace("]","")
               .replace("'","")
               .replace(" ",""))

    for suffix, niter in (('clarkclean1000', 1000), ('clarkclean10000', 10000), ):


        imagename = 'OrionFullField.{3}.spw{0}.{2}.{1}'.format(spwname, suffix, linename, band)
        if not os.path.exists("{0}.image.pbcor.fits".format(imagename)):
            print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
            assert ',' not in imagename
            assert '[' not in imagename
            tclean(vis=mses[band],
                   imagename=imagename,
                   field='Orion_BNKL_source_I',
                   spw=spw,
                   gridder='standard',
                   specmode='cube',
                   start='{0}km/s'.format(start_vel.to(u.km/u.s).value),
                   nchan=nchan,
                   restfreq='{0}Hz'.format(frequency),
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
                   chanchunks=8,
                  )
            makefits(imagename)

import numpy as np
import paths
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
import astropy.visualization
from astropy.convolution import convolve, Gaussian2DKernel
from mpl_plot_templates import asinh_norm
import matplotlib
from collections import defaultdict
import warnings
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

psf_center = coordinates.SkyCoord('5:35:14.511 -5:22:30.561',
                                  unit=(u.h, u.deg),
                                  frame='icrs')

zoomregions = {'SourceI':
               {'bottomleft': coordinates.SkyCoord("5:35:14.532",
                                                   "-5:22:30.810",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.502",
                                                 "-5:22:30.410",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.25,0.9],
                'loc': 2,
                'l1':3,
                'l2':1,
                'min': -0.0001,
                'max': 0.0002,
                'zoom': 10,
               },
               'BN':
               {'bottomleft': coordinates.SkyCoord("5:35:14.120",
                                                   "-5:22:22.820",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.092",
                                                 "-5:22:22.46",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.725,0.85],
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -0.0001,
                'max': 0.0002,
                'zoom': 10,
               },
               'WMJ053514.797-052230.557':
               {'bottomleft': coordinates.SkyCoord("5:35:14.816",
                                                   "-5:22:30.84",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.790",
                                                 "-5:22:30.49",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.25,0.3],
                'loc': 2,
                'l1':2,
                'l2':4,
                'min': -0.0001,
                'max': 0.0002,
                'zoom': 10,
               },
               'SourceN':
               {'bottomleft': coordinates.SkyCoord("5:35:14.372",
                                                   "-5:22:32.94",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.345",
                                                 "-5:22:32.50",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.75,0.3],
                'loc': 2,
                'l1':1,
                'l2':3,
                'min': -0.0001,
                'max': 0.0002,
                'zoom': 10,
               },
               'IRC6E':
               {'bottomleft': coordinates.SkyCoord("5:35:14.266",
                                                   "-5:22:28.186",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.249",
                                                 "-5:22:27.926",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.75,0.65],
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -0.001,
                'max': 0.01,
                'zoom': 10,
               },
               'IRC2C':
               {'bottomleft': coordinates.SkyCoord("5:35:14.405",
                                                   "-5:22:30.51",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.394",
                                                 "-5:22:30.32",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.75,0.50],
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -0.001,
                'max': 0.01,
                'zoom': 10,
               },
               'binary':
               {'bottomleft': coordinates.SkyCoord("5:35:14.435",
                                                   "-5:22:28.55",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.403",
                                                 "-5:22:28.28",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.5,0.85],
                'loc': 2,
                'l1':3,
                'l2':4,
                'min': -0.001,
                'max': 0.01,
                'zoom': 10,
               },
               'hotcoredisk':
               {'bottomleft': coordinates.SkyCoord("5:35:14.585",
                                                   "-5:22:31.44",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("5:35:14.564",
                                                 "-5:22:31.13",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'SourceI',
                'bbox':[0.5,0.25],
                'loc': 2,
                'l1':1,
                'l2':2,
                'min': -0.001,
                'max': 0.01,
                'zoom': 10,
               },
              }

def inset_overlays(fn, zoomregions, fignum=1,
                   psffn=None,
                   vmin=-0.001, vmax=0.01,
                   bottomleft=coordinates.SkyCoord('5:35:15.236 -5:22:39.85', unit=(u.h, u.deg), frame='icrs'),
                   topright=coordinates.SkyCoord('5:35:13.686 -5:22:19.12', unit=(u.h, u.deg), frame='icrs'),
                   tick_fontsize=pl.rcParams['axes.labelsize']):

    hdu = fits.open(fn)[0]
    print(fn)

    mywcs = wcs.WCS(hdu.header).celestial

    figure = pl.figure(fignum)
    figure.clf()
    ax = figure.add_axes([0.15, 0.1, 0.8, 0.8], projection=mywcs)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    ra.set_axislabel("RA (ICRS)", fontsize=pl.rcParams['axes.labelsize'])
    dec.set_axislabel("Dec (ICRS)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
    ra.ticklabels.set_fontsize(tick_fontsize)
    ra.set_ticks(exclude_overlapping=True)
    dec.ticklabels.set_fontsize(tick_fontsize)
    dec.set_ticks(exclude_overlapping=True)


    im = ax.imshow(hdu.data.squeeze()*1e3,
                   transform=ax.get_transform(mywcs),
                   vmin=vmin*1e3, vmax=vmax*1e3, cmap=pl.cm.gray_r,
                   interpolation='nearest',
                   origin='lower', norm=asinh_norm.AsinhNorm())

    (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg,
                                             bottomleft.dec.deg]],0)[0],
                       mywcs.wcs_world2pix([[topright.ra.deg,
                                             topright.dec.deg]],0)[0]
                      )

    # we'll want this later
    #make_scalebar(ax, scalebarpos,
    #              length=(0.5*u.pc / distance).to(u.arcsec,
    #                                              u.dimensionless_angles()),
    #              color='k',
    #              label='0.5 pc',
    #              text_offset=1.0*u.arcsec,
    #             )


    ax.set_aspect(1)
    ax.axis([x1,x2,y1,y2])


    for zoomregion in zoomregions:

        ZR = zoomregions[zoomregion]

        parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

        bl, tr = ZR['bottomleft'],ZR['topright'],
        (zx1,zy1),(zx2,zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                     bl.dec.deg]],0)[0],
                               mywcs.wcs_world2pix([[tr.ra.deg,
                                                     tr.dec.deg]],0)[0]
                              )
        print(zoomregion,zx1,zy1,zx2,zy2)

        inset_data = hdu.data.squeeze()[int(zy1):int(zy2), int(zx1):int(zx2)]
        #inset_data = hdu.data.squeeze()
        inset_wcs = mywcs.celestial[int(zy1):int(zy2), int(zx1):int(zx2)]
        #inset_wcs = mywcs

        axins = zoomed_inset_axes(parent_ax, zoom=ZR['zoom'], loc=ZR['loc'],
                                  bbox_to_anchor=ZR['bbox'],
                                  bbox_transform=figure.transFigure,
                                  axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                  axes_kwargs=dict(wcs=inset_wcs))
        ZR['axins'] = axins
        imz = axins.imshow(inset_data,
                           #transform=parent_ax.get_transform(inset_wcs),
                           extent=[int(zx1), int(zx2), int(zy1), int(zy2)],
                           vmin=ZR['min'], vmax=ZR['max'], cmap=pl.cm.gray_r,
                           interpolation='nearest',
                           origin='lower', norm=asinh_norm.AsinhNorm())


        ax.axis([x1,x2,y1,y2])
        #axins.axis([zx1,zx2,zy1,zy2])
        #print(axins.axis())

        axins.set_xticklabels([])
        axins.set_yticklabels([])

        lon = axins.coords['ra']
        lat = axins.coords['dec']
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(parent_axes=parent_ax, inset_axes=axins,
                   loc1=ZR['l1'], loc2=ZR['l2'], fc="none", ec="0.5",
                   lw=0.5)


        figure.canvas.draw()
        assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

    cax = figure.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                           ax.bbox._bbox.y1-ax.bbox._bbox.y0])
    cb = figure.colorbar(mappable=im, cax=cax)
    #print("1. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_label("$S_{1 mm}$ [mJy beam$^{-1}$]")
    #print("2. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.formatter.format = "%3.1f"
    #print("3. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_ticks(cb.formatter.locs)
    #print("4. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.set_ticklabels(["{0:3.1f}".format(float(x)) for x in cb.formatter.locs])
    #print("5. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))
    cb.ax.set_yticklabels(["{0:3.1f}".format(float(x.get_text())) for x in cb.ax.get_yticklabels()])
    #print("6. cb labels: {0}".format([x.get_text() for x in cb.ax.get_yticklabels()]))



    if psffn is not None:
        psf = fits.open(psffn)
        psfwcs = wcs.WCS(psf[0].header)
        cx,cy = psfwcs.celestial.wcs_world2pix(psf_center.ra.deg, psf_center.dec.deg, 0)
        cx = int(cx)
        cy = int(cy)
        zy1 = cy-50
        zy2 = cy+50
        zx1 = cx-50
        zx2 = cx+50

        inset_wcs = psfwcs.celestial[zy1:zy2, zx1:zx2]
        inset_data = psf[0].data[cy-50:cy+50, cx-50:cx+50]

        axins = zoomed_inset_axes(parent_ax, zoom=10, loc=2,
                                  bbox_to_anchor=(0.05,0.25),
                                  bbox_transform=figure.transFigure,
                                  axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                  axes_kwargs=dict(wcs=inset_wcs),
                                 )
        imz = axins.imshow(inset_data,
                           extent=[int(zx1), int(zx2), int(zy1), int(zy2)],
                           vmin=0, vmax=1, cmap=pl.cm.gray_r,
                           interpolation='nearest',
                           origin='lower', norm=asinh_norm.AsinhNorm())
        axins.contour(np.linspace(zx1, zx2, inset_data.shape[1]),
                      np.linspace(zy1, zy2, inset_data.shape[0]),
                      inset_data,
                      levels=[0.05, 0.1, 0.2, 0.3],
                      linewidths=[0.3]*10,
                      alpha=0.75,
                      #colors=['r']*10,
                     )
        axins.set_xticks([])
        axins.set_yticks([])
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        lon = axins.coords['ra']
        lat = axins.coords['dec']
        lon.set_ticklabel_visible(False)
        lat.set_ticklabel_visible(False)
        lon.set_ticks_visible(False)
        lat.set_ticks_visible(False)


    return figure


if __name__ == "__main__":

    zoomregions['SourceI']['min'] = -0.001
    zoomregions['SourceI']['max'] = 0.02
    zoomregions['BN']['min'] = -0.001
    zoomregions['BN']['max'] = 0.1
    zoomregions['WMJ053514.797-052230.557']['min'] = -0.001
    zoomregions['WMJ053514.797-052230.557']['max'] = 0.01
    zoomregions['SourceN']['min'] = -0.001
    zoomregions['SourceN']['max'] = 0.005
    zoomregions['IRC6E']['min'] = -0.001
    zoomregions['IRC6E']['max'] = 0.006
    zoomregions['IRC2C']['min'] = -0.001
    zoomregions['IRC2C']['max'] = 0.007
    zoomregions['binary']['min'] = -0.0005
    zoomregions['binary']['max'] = 0.002
    zoomregions['hotcoredisk']['min'] = -0.001
    zoomregions['hotcoredisk']['max'] = 0.004

    for fn in (#"Orion_SourceI_B6_continuum_r-2_longbaselines.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2_automultithresh_1mJy.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r0.5.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean5mJy.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean4mJy.selfcal.phase0.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean3mJy.selfcal.phase1.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean2mJy.selfcal.phase2.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean1mJy.selfcal.phase3.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean5mJy.automultithresh.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean4mJy.automultithresh.selfcal.phase0.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean3mJy.automultithresh.selfcal.phase1.image.tt0.pbcor.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean2mJy.automultithresh.selfcal.phase2.image.tt0.pbcor.fits",
               "Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.ampphase5.deepmask.image.tt0.pbcor.fits",
               "Orion_SourceI_B7_continuum_r-2.clean0.1mJy.500klplus.deepmask.image.tt0.pbcor.fits",
               "Orion_SourceI_B3_continuum_r-2.clean0.25mJy.image.tt0.pbcor.fits",
              ):
        figure = inset_overlays(fn, zoomregions=zoomregions,
                                psffn=fn.replace("image.tt0.pbcor","psf.tt0"),
                                vmin=-0.0005, vmax=0.005)
        figure.savefig(fn.replace(".fits","_inset.pdf"), bbox_inches='tight', dpi=300)




    zoomregions['SourceI']['min'] = -0.001
    zoomregions['SourceI']['max'] = 0.004
    zoomregions['BN']['min'] = -0.001
    zoomregions['BN']['max'] = 0.004
    zoomregions['WMJ053514.797-052230.557']['min'] = -0.001
    zoomregions['WMJ053514.797-052230.557']['max'] = 0.004
    zoomregions['SourceN']['min'] = -0.001
    zoomregions['SourceN']['max'] = 0.004
    zoomregions['IRC6E']['min'] = -0.001
    zoomregions['IRC6E']['max'] = 0.004
    zoomregions['IRC2C']['min'] = -0.001
    zoomregions['IRC2C']['max'] = 0.004
    zoomregions['binary']['min'] = -0.001
    zoomregions['binary']['max'] = 0.001
    zoomregions['hotcoredisk']['min'] = -0.001
    zoomregions['hotcoredisk']['max'] = 0.001

    for fn in (#"Orion_SourceI_B6_continuum_r-2_longbaselines.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2_automultithresh_1mJy.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r0.5.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean5mJy.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean4mJy.selfcal.phase0.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean3mJy.selfcal.phase1.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean2mJy.selfcal.phase2.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean1mJy.selfcal.phase3.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean0.5mJy.selfcal.phase4.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean5mJy.automultithresh.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean4mJy.automultithresh.selfcal.phase0.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean3mJy.automultithresh.selfcal.phase1.residual.tt0.fits",
               #"Orion_SourceI_B6_continuum_r-2.clean2mJy.automultithresh.selfcal.phase2.residual.tt0.fits",
               "Orion_SourceI_B3_continuum_r-2.clean0.25mJy.residual.tt0.fits",
               "Orion_SourceI_B6_continuum_r-2.clean0.1mJy.selfcal.ampphase5.deepmask.residual.tt0.fits",
              ):
        figure = inset_overlays(fn, zoomregions=zoomregions,
                                psffn=fn.replace("residual","psf"),
                                vmin=-0.0001, vmax=0.001)
        figure.savefig(fn.replace(".fits","_inset.pdf"), bbox_inches='tight',
                       dpi=300)

# Orion_SourceI_B3_continuum_r0.5_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r0.5.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r-2_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r2_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r-2.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r2.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r-2.mask2.5mJy.clean0.5mJy.image.tt0.pbcor.fits
# Orion_SourceI_B3_continuum_r-2.mask2mJy.clean1mJy.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r0.5_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r0.5.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2_automultithresh_1mJy.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r2_dirty.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r2.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2_longbaselines.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.selfcal.phase0.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.selfcal.phase1.image.tt0.pbcor.fits
# Orion_SourceI_B6_continuum_r-2.mask5mJy.clean4mJy.selfcal.phase2.image.tt0.pbcor.fits

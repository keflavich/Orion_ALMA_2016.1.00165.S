import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import coordinates, units as u
from astropy import visualization
from astropy.visualization import simple_norm
import regions

import pylab as pl

from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, BboxConnector
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes

from astropy import units as u
from astropy import coordinates
from astropy.stats import mad_std
from astropy.visualization import simple_norm
import astropy.visualization
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import radio_beam

import paths

distance = 400*u.pc

fontsize = 16

def make_scalebar(ax, left_side, length, color='w', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform('fk5'),
                   )
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform('fk5'),
                  color=color,
                  fontsize=fontsize,
                 )
    ax.axis(axlims)
    return lines,txt

def hide_ticks(ax):
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ra.set_ticks_visible(False)
    dec.set_ticks_visible(False)
    ra.set_axislabel('')
    dec.set_axislabel('')
    ra.ticklabels.set_visible(False)
    dec.ticklabels.set_visible(False)

def mark_inset_otherdata(axins, parent_ax, bl, tr, loc1, loc2, edgecolor='b', **markinkwargs):
    blt = bl.transform_to(parent_ax.wcs.wcs.radesys.lower())
    trt = tr.transform_to(parent_ax.wcs.wcs.radesys.lower())
    (rx1,ry1),(rx2,ry2) = (parent_ax.wcs.wcs_world2pix([[blt.ra.deg,
                                                         blt.dec.deg]],0)[0],
                           parent_ax.wcs.wcs_world2pix([[trt.ra.deg,
                                                         trt.dec.deg]],0)[0]
                          )
    bbox = Bbox(np.array([(rx1,ry1),(rx2,ry2)]))
    rect = TransformedBbox(bbox, parent_ax.transData)

    markinkwargs = dict(fc='none', ec=edgecolor, **markinkwargs)

    pp = BboxPatch(rect, fill=False, **markinkwargs)
    parent_ax.add_patch(pp)

    p1 = BboxConnector(axins.bbox, rect, loc1=loc1, **markinkwargs)
    axins.add_patch(p1)
    p1.set_clip_on(False)
    p2 = BboxConnector(axins.bbox, rect, loc1=loc2, **markinkwargs)
    axins.add_patch(p2)
    p2.set_clip_on(False)

    return bbox, rect, p1, p2


for largescale, largescalename, stretch in (('/Users/adam/work/bolocam/orion/orionfruitawf2.fits', 'continuum', 'asinh'),
                                   ('/Users/adam/work/orion/alma/ISF_ALMA+IRAM30m_N2H+10_mom0.fits', 'n2hp', 'asinh'),
                                   ('/Users/adam/work/orion/oriona_isf_artemis_coldens.fit', 'artemisNH2', 'asinh'),
                                   ('/Users/adam/work/orion/oriona_isf_artemis_temp.fit', 'artemisT', 'linear'),
                                  ):


    fh1 = fits.open(largescale)
    fh2 = fits.open('/Users/adam/work/orion/alma/FITS/Orion_NW_12m_7m_merge_continuum.fits')
    fh3 = fits.open('/Users/adam/Dropbox/Orion_ALMA_LB/FITS/Orion_SourceI_B7_continuum_r2.clean1mJy.allbaselines.deepmask.image.tt0.pbcor.fits')
    fh4 = fits.open('/Users/adam/work/orion/alma_lb/FITS/Orion_SourceI_B7_continuum_r-2.clean0.1mJy.500klplus.deepmask.image.tt0.pbcor.fits')
    w1 = wcs.WCS(fh1[0].header).celestial
    w2 = wcs.WCS(fh2[0].header)
    w3 = wcs.WCS(fh3[0].header)
    w4 = wcs.WCS(fh4[0].header)

    zoom1,zoom2,zoom3,zoom4,zoom5 = regions.Regions.parse("""fk5
    box(5:35:13.4058,-5:22:28.223,448.906",590.665",359.99999)
    box(5:35:14.2200,-5:22:25.065,45",53.734",359.99999)
    box(5:35:14.4600,-5:22:31.506,10.0",13.034",359.99999)
    box(5:35:14.5772,-5:22:31.319,0.443",0.40",0)
    box(5:35:14.4218,-5:22:28.451,0.443",0.443",3.2657584e-07)
                                                    """, format='ds9')
    # SrcI box(5:35:14.5187,-5:22:30.653,0.443",0.443",0)


    msk1 = zoom1.to_pixel(w1).to_mask()
    w1z = w1[msk1.get_overlap_slices(fh1[0].data.squeeze().shape)[0]]

    fig = pl.figure(1, figsize=(12,10))
    fig.clf()

    ax = fig.add_subplot(position=(0.1,0.1,0.4,0.8), projection=w1z)

    data1 = msk1.multiply(fh1[0].data.squeeze())
    ax.imshow(data1,  cmap='magma',
                  norm=simple_norm(data1, stretch=stretch)
    )

    blx, urx, bly, ury = ax.axis([zoom1.to_pixel(w1z).center.x-zoom1.to_pixel(w1z).width/2,
             zoom1.to_pixel(w1z).center.x+zoom1.to_pixel(w1z).width/2,
             zoom1.to_pixel(w1z).center.y-zoom1.to_pixel(w1z).height/2,
             zoom1.to_pixel(w1z).center.y+zoom1.to_pixel(w1z).height/2,
            ])

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    dec.set_major_formatter('dd:mm:ss.s')
    ra.ticklabels.set_fontsize(16)
    dec.ticklabels.set_fontsize(16)

    ax.set_xlabel('Right Ascension', fontsize=16)
    ax.set_ylabel('Declination', fontsize=16)


    scalebar_length = 0.1*u.pc
    length = (scalebar_length / distance).to(u.arcsec, u.dimensionless_angles())
    bls = w1z.pixel_to_world(blx+zoom1.to_pixel(w1z).width*5/6, bly+zoom1.to_pixel(w1z).height/10)
    left_side = coordinates.SkyCoord(bls, frame=w1z.wcs.radesys.lower(), unit=(u.deg,u.deg))
    make_scalebar(ax, left_side, length, color='w', linestyle='-', label=f'{scalebar_length:0.1f}',
                  text_offset=0.1*u.arcsec, fontsize=fontsize)


    msk2 = zoom2.to_pixel(w2).to_mask()
    w2z = w2[msk2.get_overlap_slices(fh2[0].data.shape)[0]]
    axins2 = inset_axes(ax,
                       loc=1, width=3, height=3,
                       bbox_to_anchor=(0.75,0.9),
                       bbox_transform=fig.transFigure,
                       axes_class=visualization.wcsaxes.core.WCSAxes,
                       axes_kwargs=dict(wcs=w2z))

    data2 = msk2.multiply(fh2[0].data)

    axins2.imshow(data2,  cmap='inferno',
                  norm=simple_norm(data2, stretch='asinh')
                 )

    blx,urx,bly,ury = axins2.axis([zoom2.to_pixel(w2z).center.x-zoom2.to_pixel(w2z).width/2,
                                   zoom2.to_pixel(w2z).center.x+zoom2.to_pixel(w2z).width/2,
                                   zoom2.to_pixel(w2z).center.y-zoom2.to_pixel(w2z).height/2,
                                   zoom2.to_pixel(w2z).center.y+zoom2.to_pixel(w2z).height/2,
                                   ])
    if urx < blx:
        urx, blx = blx, urx
    bl = w2z.pixel_to_world(blx, bly)
    tr = w2z.pixel_to_world(urx, ury)
    #ax.plot(bl.ra, bl.dec, transform=ax.get_transform('fk5'), color='r', marker='x', linestyle='none')
    #ax.plot(tr.ra, tr.dec, transform=ax.get_transform('fk5'), color='c', marker='o', linestyle='none')
    mark_inset_otherdata(axins2, ax, bl, tr, 2, 4,)
    hide_ticks(axins2)


    scalebar_length = 5000*u.au
    length = (scalebar_length / distance).to(u.arcsec, u.dimensionless_angles())
    bls = w2z.pixel_to_world(blx+zoom2.to_pixel(w2z).width/4, bly)
    left_side = coordinates.SkyCoord(bls, frame=w2z.wcs.radesys.lower(), unit=(u.deg,u.deg))
    make_scalebar(axins2, left_side, length, color='w', linestyle='-', label=f'{scalebar_length:0.1f}',
                  text_offset=0.1*u.arcsec, fontsize=fontsize)
    #axins2.text(0.05, 0.9, 'Inset 1', transform=axins2.transAxes, color='w', fontsize=fontsize)



    msk3 = zoom3.to_pixel(w3).to_mask()
    w3z = w3[msk3.get_overlap_slices(fh3[0].data.shape)[0]]
    axins3 = inset_axes(ax,
                       loc=4, width=3, height=3,
                       bbox_to_anchor=(0.75,0.1),
                       bbox_transform=fig.transFigure,
                       axes_class=visualization.wcsaxes.core.WCSAxes,
                       axes_kwargs=dict(wcs=w3z))

    data3 = msk3.multiply(fh3[0].data)

    axins3.imshow(data3,  cmap='inferno',
                  norm=simple_norm(data3, stretch='asinh', max_percent=99.99)
                 )

    blx,urx,bly,ury = axins3.axis([zoom3.to_pixel(w3z).center.x-zoom3.to_pixel(w3z).width/2,
                                   zoom3.to_pixel(w3z).center.x+zoom3.to_pixel(w3z).width/2,
                                   zoom3.to_pixel(w3z).center.y-zoom3.to_pixel(w3z).height/2,
                                   zoom3.to_pixel(w3z).center.y+zoom3.to_pixel(w3z).height/2,
                                   ])
    if urx < blx:
        urx, blx = blx, urx
    bl = w3z.pixel_to_world(blx, bly)
    tr = w3z.pixel_to_world(urx, ury)
    #ax.plot(bl.ra, bl.dec, transform=ax.get_transform('fk5'), color='r', marker='x', linestyle='none')
    #ax.plot(tr.ra, tr.dec, transform=ax.get_transform('fk5'), color='c', marker='o', linestyle='none')
    mark_inset_otherdata(axins3, axins2, bl, tr, 1, 2,)
    hide_ticks(axins3)


    scalebar_length = 500*u.au
    length = (scalebar_length / distance).to(u.arcsec, u.dimensionless_angles())
    bls = w3z.pixel_to_world(blx+zoom3.to_pixel(w3z).width/4, bly)
    left_side = coordinates.SkyCoord(bls, frame=w3z.wcs.radesys.lower(), unit=(u.deg,u.deg))
    make_scalebar(axins3, left_side, length, color='w', linestyle='-', label=f'{scalebar_length:0.1f}',
                  text_offset=0.1*u.arcsec, fontsize=fontsize)
    #axins3.text(0.05, 0.9, 'Inset 2', transform=axins3.transAxes, color='w', fontsize=fontsize)





    msk4 = zoom4.to_pixel(w4).to_mask()
    w4z = w4[msk4.get_overlap_slices(fh4[0].data.shape)[0]]
    axins4 = inset_axes(ax,
                       loc=1, width=3, height=3,
                       bbox_to_anchor=(1.0,0.9),
                       bbox_transform=fig.transFigure,
                       axes_class=visualization.wcsaxes.core.WCSAxes,
                       axes_kwargs=dict(wcs=w4z))

    data4 = msk4.multiply(fh4[0].data)

    axins4.imshow(data4,  cmap='inferno',
                  norm=simple_norm(data4, stretch='asinh', max_percent=99.99, min_percent=1)
                 )

    blx,urx,bly,ury = axins4.axis([zoom4.to_pixel(w4z).center.x-zoom4.to_pixel(w4z).width/2,
                                   zoom4.to_pixel(w4z).center.x+zoom4.to_pixel(w4z).width/2,
                                   zoom4.to_pixel(w4z).center.y-zoom4.to_pixel(w4z).height/2,
                                   zoom4.to_pixel(w4z).center.y+zoom4.to_pixel(w4z).height/2,
                                   ])
    if urx < blx:
        urx, blx = blx, urx
    bl = w4z.pixel_to_world(blx, bly)
    tr = w4z.pixel_to_world(urx, ury)
    #ax.plot(bl.ra, bl.dec, transform=ax.get_transform('fk5'), color='r', marker='x', linestyle='none')
    #ax.plot(tr.ra, tr.dec, transform=ax.get_transform('fk5'), color='c', marker='o', linestyle='none')
    _,_,p1,p2 = mark_inset_otherdata(axins4, axins3, bl, tr, 4, 2, edgecolor='g')
    hide_ticks(axins4)


    scalebar_length = 50*u.au
    length = (scalebar_length / distance).to(u.arcsec, u.dimensionless_angles())
    bls = w4z.pixel_to_world(blx+zoom4.to_pixel(w4z).width/4, bly)
    left_side = coordinates.SkyCoord(bls, frame=w4z.wcs.radesys.lower(), unit=(u.deg,u.deg))
    make_scalebar(axins4, left_side, length, color='w', linestyle='-', label=f'{scalebar_length:0.1f}',
                  text_offset=0.01*u.arcsec, fontsize=fontsize)
    #axins4.set_zorder(-5)
    #axins4.text(0.05, 0.9, 'Inset 3', transform=axins4.transAxes, color='w', fontsize=fontsize, text_offset=0.001*u.arcsec)




    msk5 = zoom5.to_pixel(w4).to_mask()
    w5z = w4[msk5.get_overlap_slices(fh4[0].data.shape)[0]]
    axins5 = inset_axes(axins3,
                       loc=1, width=3, height=3,
                       bbox_to_anchor=(1.0,0.5),
                       bbox_transform=fig.transFigure,
                       axes_class=visualization.wcsaxes.core.WCSAxes,
                       axes_kwargs=dict(wcs=w5z))

    data5 = msk5.multiply(fh4[0].data)

    axins5.imshow(data5,  cmap='inferno',
                  norm=simple_norm(data5, stretch='asinh', max_percent=99.99, min_percent=1)
                 )

    blx,urx,bly,ury = axins5.axis([zoom5.to_pixel(w5z).center.x-zoom5.to_pixel(w5z).width/2,
                                   zoom5.to_pixel(w5z).center.x+zoom5.to_pixel(w5z).width/2,
                                   zoom5.to_pixel(w5z).center.y-zoom5.to_pixel(w5z).height/2,
                                   zoom5.to_pixel(w5z).center.y+zoom5.to_pixel(w5z).height/2,
                                   ])
    if urx < blx:
        urx, blx = blx, urx
    bl = w5z.pixel_to_world(blx, bly)
    tr = w5z.pixel_to_world(urx, ury)
    #ax.plot(bl.ra, bl.dec, transform=ax.get_transform('fk5'), color='r', marker='x', linestyle='none')
    #ax.plot(tr.ra, tr.dec, transform=ax.get_transform('fk5'), color='c', marker='o', linestyle='none')
    hide_ticks(axins5)
    _,_,p1,p2 = mark_inset_otherdata(axins5, axins3, bl, tr, 2, 3, edgecolor='r')


    scalebar_length = 50*u.au
    length = (scalebar_length / distance).to(u.arcsec, u.dimensionless_angles())
    bls = w5z.pixel_to_world(blx+zoom5.to_pixel(w5z).width/5, bly)
    left_side = coordinates.SkyCoord(bls, frame=w5z.wcs.radesys.lower(), unit=(u.deg,u.deg))
    make_scalebar(axins5, left_side, length, color='w', linestyle='-', label=f'{scalebar_length:0.1f}',
                  text_offset=0.01*u.arcsec, fontsize=fontsize)

    pl.savefig(paths.fpath(f"HotCoreDisk_ZoomIn_{largescalename}.png"), bbox_inches='tight')
    pl.savefig(paths.fpath(f"HotCoreDisk_ZoomIn_{largescalename}.pdf"), bbox_inches='tight')

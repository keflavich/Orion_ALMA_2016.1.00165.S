import numpy as np
from astropy import stats
from scipy.ndimage.measurements import label

def mad_std_nan(*args, **kwargs):
    return stats.mad_std(*args, ignore_nan=True, **kwargs)

def find_continuum_channels(spectrum, stdfunc=mad_std_nan):
    """
    Give a spectrum, run astropy's sigma-clipping algorithm to identify the
    continuum channels.  Return a label array that will be nonzero wherever
    the channels are either lines or NaNs and zero wherever the data are
    most likely continuum.
    """

    masked_data = stats.sigma_clip(spectrum, stdfunc=stdfunc)

    continua = (~masked_data.mask)

    # continuum labels (not generally needed)
    # labels, nlabels = (label(continua))
    # Llabels = line labels
    Llabels, nLlabels = (label(~continua))

    return Llabels, nLlabels, masked_data

def casaify_labels(labels, nlabels, spw, xarr=None, min_flag_length=0):
    """
    Given a set of labels (assignments of an array to groups), print out the
    CASA version of both the per-channel and the by-frequency (if xarr is
    provided) masks for flagging those data out
    """

    bands = ''
    fbands = ''

    for ii in range(1, nlabels+1):

        # select the stuff we want to flag out
        x = np.where(labels == ii)
        if len(x[0]) > min_flag_length:
            bands = bands + ",{spw}:{0}~{1}".format(x[0][0], x[0][-1], spw=spw)
            if xarr is not None:
                fbands = (fbands +
                          ",{spw}:{0}~{1}GHz".format(xarr[x[0][0]].value,
                                                     xarr[x[0][-1]].value,
                                                     spw=spw)
                         )

    return bands.lstrip(","), fbands.lstrip(",")

def casamask_to_boolmask(casamask, size=None):

    masks = {}

    for spwchunk in casamask.split(","):
        colonsplit = spwchunk.split(":")
        if len(colonsplit) > 1:
            spw = int(colonsplit[0])
            chans = colonsplit[1]
        else:
            chans = colonsplit[0]

        ranges = [list(map(int, xx.split("~"))) for xx in chans.split(";")]
        maxval = np.max(ranges)
        if size is None:
            size = maxval
        boolarr = np.zeros(size, dtype='bool')
        for low,high in ranges:
            boolarr[low:high] = True

        masks[spw] = boolarr

    return masks


if __name__ == "__main__":

    # use the spectra from average_spectra to identify line / continuum channels
    
    from spectral_cube.lower_dimensional_structures import OneDSpectrum
    from astropy.io import fits
    import paths
    import pylab as pl

    qa2_flagchannels='0:80~110;155~195;235~260;325~355;420~440;480~495;585~610;695~720;770~790;865~885;1215~1230;1310~1495;1595~1615,1:0~170;315~450;585~600;670~715;845~890;960~985;1125~1875,2:35~90;305~345;470~490;530~605;660~780;865~1150;1260~1285;1450~1565;1685~1705;1785~1870,3:80~125;235~660;985~995;1045~1075;1235~1245;1585~1670;1765~1830'
    qa2_masks = casamask_to_boolmask(qa2_flagchannels, size=1920)

    pl.figure(1).clf()

    ftemplate = 'OrionSourceI_only.B6.robust0.5.spw{0}.maskedclarkclean10000_medsub_avg.fits'

    fullbands = ''
    for spw in (0,1,2,3):

        fn = paths.dpath('spectra/'+ftemplate.format(spw))
        avgspec = OneDSpectrum.from_hdu(fits.open(fn))

        labels, nlabels, masked_data = find_continuum_channels(avgspec.value)
        bands, fbands = casaify_labels(labels, nlabels, spw=spw, min_flag_length=1)

        fullbands += bands
        print(bands)

        ax = pl.subplot(4,1,spw+1)
        ax.plot(avgspec.spectral_axis, avgspec.value, color='k')
        ax.plot(avgspec.spectral_axis, masked_data, color='r', linewidth=3,
                alpha=0.5, zorder=10)
        ax.plot(avgspec.spectral_axis, np.ma.masked_where(qa2_masks[spw],
                                                          avgspec.value),
                color='b', linewidth=5,
                alpha=0.35, zorder=15)

    with open(paths.redpath('linechannels_b6'), 'w') as fh:
        fh.write(fullbands)

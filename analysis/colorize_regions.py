import regions
import paths
import pylab as pl

for linename,(vmin,vmax),limits in (('Unknown_4', (-15, 27), (-0.1, 0.1, -0.12, 0.12)),
                                    ('SiOv=1_5-4', (-23, 30), (-0.2, 0.2, -0.2, 0.2)),
                                   ):

    regs = regions.read_ds9(paths.rpath('velo_centroid_guesses_{linename}.reg').format(linename=linename))

    cmap = pl.cm.Spectral_r
    cmap = pl.cm.spectral
    norm = pl.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    with open(paths.rpath('velo_centroid_guesses_{linename}.reg').format(linename=linename), 'w') as fh:
        fh.write('image\n')

        for reg in regs:
            text = reg.meta['text'].strip('{}')
            velo = float(text)
            color = pl.matplotlib.colors.rgb2hex(cmap(norm(velo)))
            fh.write("point({0}, {1}) # point=x color={color} text={{{velo}}}\n"
                     .format(reg.center.x, reg.center.y, color=color,
                             velo=velo))

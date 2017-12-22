import regions
import paths
import pylab as pl

regs = regions.read_ds9(paths.rpath('velo_centroid_guesses_Unknown_4.reg'))

cmap = pl.cm.Spectral_r
cmap = pl.cm.spectral
vmin, vmax = -15, 27
norm = pl.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

with open(paths.rpath('velo_centroid_guesses_Unknown_4.reg'), 'w') as fh:
    fh.write('image\n')

    for reg in regs:
        text = reg.meta['text'].strip('{}')
        velo = float(text)
        color = pl.matplotlib.colors.rgb2hex(cmap(norm(velo)))
        fh.write("point({0}, {1}) # point=x color={color} text={{{velo}}}\n"
                 .format(reg.center.x, reg.center.y, color=color,
                         velo=velo))

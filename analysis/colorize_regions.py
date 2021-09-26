import regions
import paths
import pylab as pl

def colorize_region(fn, vlims=None, ):

    regs = regions.Regions.read(fn)

    if vlims is not None:
        vmin,vmax = vlims
    else:
        vels = [float(reg.meta['text'].strip("{}"))
                for reg in regs]
        vmin, vmax = min(vels), max(vels)

    cmap = pl.cm.Spectral_r
    cmap = pl.cm.spectral
    norm = pl.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

    for reg in regs:
        text = reg.meta['text'].strip('{}')
        velo = float(text)
        color = pl.matplotlib.colors.rgb2hex(cmap(norm(velo)))
        reg.meta['color'] = color

    return regs

if __name__ == "__main__":
    for linename,(vmin,vmax),limits in (('Unknown_4', (-15, 27), (-0.1, 0.1, -0.12, 0.12)),
                                        ('Unknown_1', (-15, 27), (-0.1, 0.1, -0.12, 0.12)),
                                        ('SiOv=1_5-4', (-30, 45), (-0.2, 0.2, -0.2, 0.2)),
                                        ('H2Ov2=1_5(5,0)-6(4,3)', (-28, 38), (-0.2, 0.2, -0.2, 0.2)),
                                       ):

        regname = paths.rpath('velo_centroid_guesses_{linename}.reg').format(linename=linename)

        regs = colorize_region

        with open(paths.rpath('velo_centroid_guesses_{linename}.reg').format(linename=linename), 'w') as fh:
            fh.write('image\n')

            for reg in regs:
                text = reg.meta['text'].strip('{}')
                velo = float(text)
                color = reg.meta['color'].strip('{}')
                fh.write("point({0}, {1}) # point=x color={color} text={{{velo}}}\n"
                         .format(reg.center.x, reg.center.y, color=color,
                                 velo=velo))

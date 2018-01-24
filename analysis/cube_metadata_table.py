import paths
import os
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from spectral_cube import SpectralCube, wcs_utils, tests, Projection
import pylab as pl
from lines import disk_lines
import radio_beam


tabletext = (r"""
\begin{table*}[htp]
\centering
\caption{Line Cube Parameters}
\begin{tabular}{ccccccccc}
\label{tab:cube_metadata}
Band & SPW & Freq. Range & Robust & Beam Major & Beam Minor & Beam PA               & RMS                               & RMS\\
     &     & GHz         &        & \arcsec    & \arcsec    & $\mathrm{{}^{\circ}}$ & $\mathrm{mJy}~\mathrm{beam}^{-1}$ & K\\
\hline

XXDATAXX

\hline
\end{tabular}

\end{table*}
""")

datatext = []

for robust in (-2, 0.5):
    ftemplate = '/Volumes/external/orion/Orion{1}_only.{2}.robust{robust}.spw{0}.{suffix}_medsub.image.pbcor.fits'

    for band,suffix in (('B3', 'clarkclean10000'),
                        ('B6', 'maskedclarkclean10000')):
        for sourcename in ('SourceI',):# 'BN'):

            for spw in (0,1,2,3):
                filename = ftemplate.format(spw, sourcename, band,
                                            suffix=suffix, robust=robust)

                if os.path.exists(filename):
                    cube = SpectralCube.read(filename).mask_out_bad_beams(0.1)
                else:
                    log.exception("File {0} does not exist".format(filename))
                    continue

                fmin, fmax = cube.spectral_extrema

                noise = cube.mad_std()
                average_beam = cube.average_beams(0.1)

                row = (r"{band} & {spw} & {freq1:0.3f}-{freq2:0.3f} & "
                       r"{robust} & {bmaj:0.3f} & {bmin:0.3f} & {bpa:0.1f} & "
                       r"{rms1:0.1f} & {rms2:0.1f} \\"
                       .format(band=band,
                               bmaj=average_beam.major.to(u.arcsec).value,
                               bmin=average_beam.minor.to(u.arcsec).value,
                               bpa=average_beam.pa.to(u.deg).value,
                               rms1=(noise.to(u.mJy/u.beam)).value,
                               rms2=noise.to(u.K, average_beam.jtok_equiv((fmin+fmax)/2)).value,
                               spw=spw,
                               robust=robust,
                               freq1=fmin.to(u.GHz).value,
                               freq2=fmax.to(u.GHz).value,
                              )
                      )
                datatext.append(row)

with open(paths.texpath('cube_metadata.tex'), 'w') as texfh:

    texfh.write(tabletext.replace("XXDATAXX", "\n".join(datatext)))

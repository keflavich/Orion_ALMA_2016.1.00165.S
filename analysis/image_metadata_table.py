# make metadata table for images used in analysis

from files import b6_hires_cont, b3_hires_cont
import numpy as np
import paths
from astropy.io import fits
from astropy import units as u
from astropy import wcs
import regions
from radio_beam import Beam


noiseregion = regions.Regions.read(paths.rpath('noise_estimate_region.reg'))[0]
sourceIcircle = regions.Regions.read(paths.rpath('sourceI_enclosing_circle.reg'))[0]


tabletext = (r"""
\begin{table*}[htp]
\centering
\caption{Continuum Image Parameters}
\begin{tabular}{ccccccccc}
\label{tab:image_metadata}
Band & Robust & Beam Major & Beam Minor & Beam PA               & $T_B$/$S_\nu$      & RMS & Source I $S_{\nu,max}$ & Dynamic Range\\
     &        & \arcsec    & \arcsec    & $\mathrm{{}^{\circ}}$ & $10^3$ K Jy$^{-1}$ & $\mathrm{mJy}~\mathrm{beam}^{-1}$ & $\mathrm{mJy}~\mathrm{beam}^{-1}$ & \\
\hline

XXDATAXX

\hline
\end{tabular}

\end{table*}
""")

datatext = []

for fn,freq in ((b6_hires_cont, 224*u.GHz),
                (b3_hires_cont, 93.3*u.GHz)):

    fh = fits.open(paths.dpath(fn))
    header = fh[0].header
    data = fh[0].data
    beam = Beam.from_fits_header(header)
    ww = wcs.WCS(header)

    nrp = noiseregion.to_pixel(ww)
    mask = nrp.to_mask()
    noise = mask.cutout(data)[mask.data.astype('bool')].std()

    sourceImask = sourceIcircle.to_pixel(ww).to_mask()
    sourceIpeak = sourceImask.cutout(data)[sourceImask.data.astype('bool')].max()

    row = (r"{band} & {robust} & {bmaj:0.3f} & {bmin:0.3f} & {bpa:0.1f} &"
           r" {jtok:0.1f} &"
           r" {rms:0.3f} & {srcipeak:0.3f} & {DR:0d} \\"
           .format(band="B3" if "B3" in fn else "B6",
                   bmaj=beam.major.to(u.arcsec).value,
                   bmin=beam.minor.to(u.arcsec).value,
                   bpa=beam.pa.to(u.deg).value,
                   rms=noise*1e3,
                   robust=('-2' if 'r-2' in fn
                           else '0.5' if 'r0.5' in fn
                           else '2' if 'r2' in fn
                           else 'ERROR'),
                   srcipeak=sourceIpeak*1e3,
                   DR=int((sourceIpeak / noise)/10)*10,
                   jtok=(beam.jtok(freq).value/1e3),
                  )
          )
    datatext.append(row)



with open(paths.texpath('image_metadata.tex'), 'w') as texfh:

    texfh.write(tabletext.replace("XXDATAXX", "\n".join(datatext)))

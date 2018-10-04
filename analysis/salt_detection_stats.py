import numpy as np
import paths
from astropy.table import Table,Column
import latex_info
import re

tbl = Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')

detections = (tbl['Flag'] == '-d')

print("Counts of unambiguous detections in all isotopologues for each vibrational state:")
for vibstate in range(11):
    mask = np.array(['v={0}-{0}'.format(vibstate) in x for x in tbl['Species']])
    print("v={0} -> {1} detections".format(vibstate, (mask & detections).sum()))
ignore_vib = np.array(['v=10-10' in x for x in tbl['Species']]) | np.array(['v=9-9' in x for x in tbl['Species']])

assert not np.any(detections & ignore_vib)

tbl = tbl[~ignore_vib]

detections = np.array([x[1]=='d' for x in tbl['Flag']])
questionable = (tbl['Flag'] == '-q') | (tbl['Flag'] == 'cq')
confused = np.array([x[0]=='c' for x in tbl['Flag']])
ndetections = detections.sum()
ninband = len(tbl)

print()
print("Total of {0} unambiguous detections out of {1} lines with v<=8 in band ({2:0.1f} %)".format(ndetections, ninband, ndetections/ninband*100))
print("{0} lines were confused.  Out of the {1} unconfused lines, {2} were detected ({3:0.1f} %)"
      .format(confused.sum(), (ninband-confused.sum()), (detections & (~confused)).sum(),
              (detections & (~confused)).sum()/(ninband-confused.sum())*100))

for species in ('23Na-35Cl', '23Na-37Cl', '39K-35Cl', '41K-35Cl', '39K-37Cl', '41K-37Cl'):
    mask = np.array([species in x for x in tbl['Species']])

    ndet = (mask & detections).sum()

    print("Species {0}: {1} unambiguous detections out of {2} lines in band ({3:0.1f} %)".format(species, ndet, mask.sum(), ndet/mask.sum()*100))

print()
print("If we exclude confused lines from the count: ")
for species in ('23Na-35Cl', '23Na-37Cl', '39K-35Cl', '41K-35Cl', '39K-37Cl', '41K-37Cl'):
    mask = np.array([species in x for x in tbl['Species']])

    mask &= (~confused)
    ndet = (mask & detections).sum()

    print("Species {0}: {1} unambiguous detections out of {2} lines in band ({3:0.1f} %)".format(species, ndet, mask.sum(), ndet/mask.sum()*100))



vmax = 5
vmask = np.array([any(['v={0}-{0}'.format(v) in x for v in range(vmax+1)]) for x in tbl['Species']])
tbl = tbl[vmask]

detections = np.array([x[1]=='d' for x in tbl['Flag']])
questionable = (tbl['Flag'] == '-q') | (tbl['Flag'] == 'cq')
confused = np.array([x[0]=='c' for x in tbl['Flag']])
ndetections = detections.sum()
ninband = len(tbl)

print()
print("Total of {0} unambiguous detections out of {1} lines with v<={3} in band ({2:0.1f} %)".format(ndetections, ninband, ndetections/ninband*100, vmax))
print("{0} lines were confused.  Out of the {1} unconfused lines, {2} were detected ({3:0.1f} %)"
      .format(confused.sum(), (ninband-confused.sum()), (detections & (~confused)).sum(),
              (detections & (~confused)).sum()/(ninband-confused.sum())*100))

for species in ('23Na-35Cl', '23Na-37Cl', '39K-35Cl', '41K-35Cl', '39K-37Cl', '41K-37Cl'):
    mask = np.array([species in x for x in tbl['Species']])

    ndet = (mask & detections).sum()

    print("Species {0}: {1} unambiguous detections out of {2} lines in band ({3:0.1f} %)".format(species, ndet, mask.sum(), ndet/mask.sum()*100))



# reload the unmasked table
tbl = Table.read(paths.tpath('salts_in_band.ipac'), format='ascii.ipac')

vre = re.compile('v=([0-9]+)')
vstate = [int(vre.search(row['Species']).groups()[0]) for row in tbl]
jre = re.compile('J=([0-9]+)-([0-9]+)')
Ju,Jl = zip(*[map(int, jre.search(row['Species']).groups()) for row in tbl])
tbl.add_column(Column(name='v', data=vstate))
tbl.add_column(Column(name='J$_u$', data=Ju))
tbl.add_column(Column(name='J$_l$', data=Jl))
tbl.rename_column('E_U', 'E$_U$')
tbl.rename_column('Aij', 'A$_{ul}$')
nre = re.compile('(23|39|37|41|35)')
def expnum(x):
    return nre.sub(r'$^{\1}$', x)
species = ["".join(map(expnum, row['Species'].split('v')[0].split('-')))
           for row in tbl]
tbl.rename_column('Species', 'fullSpecies')
tbl.add_column(Column(name='Species', data=species))

latexdict = latex_info.latexdict.copy()
latexdict['preamble'] = '\centering'
latexdict['tablefoot'] = ('\n\par ')

formats = {'Frequency': lambda x: "{0:0.5f}".format(x),
           'E$_U$': lambda x: "-" if np.isnan(x) else "{0:0.1f}".format(x),
           'A$_{ul}$': lambda x: "{0:0.5f}".format(x),
          }


for band in ('B3','B6','B7.lb'):
    match = tbl['Band'] == band
    latexdict['header_start'] = '\label{{tab:all_detections_B{0}}}'.format(band[1])
    latexdict['caption'] = 'All observed lines in Band {0}'.format(band[1])
    cols = ['Species', 'v', 'J$_u$', 'J$_l$', 'E$_U$', 'A$_{ul}$', 'Frequency', 'Flag']
    print("{1}: nmatch = {0}".format(match.sum(), band))
    tbl[cols][match].write(paths.texpath2('lines_in_band{0}.tex'.format(band[1])),
                           formats=formats, latexdict=latexdict,
                           overwrite=True)

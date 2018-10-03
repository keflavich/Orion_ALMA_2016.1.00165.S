import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
from astropy import modeling
import paths
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
from astropy.table import Column
from astroquery.vizier import Vizier

def load_barton(species):
    tbl = table.Table.read(paths.salty(species+"_rotational_transitions.ipac"),
                           format='ascii.ipac')
    tbl.add_column(table.Column(name='Species', data=[species+row['QNs'] for row in tbl]))
    tbl.rename_column('Frequency', 'Freq')
    tbl = tbl[(tbl['Freq']>1) & (tbl['Freq']<400) & (tbl['E_U'] < 1e4) & (tbl['vu']<=10)]
    return tbl

kcl_offset = (1+17*u.km/u.s/constants.c).decompose()
nacl_offset = (1+3.*u.km/u.s/constants.c).decompose()
kcl_offset = 1
nacl_offset = 1



def get_offset_model(species_diff_table):
    """
    Using the overlapping species in splatalogue & Barton, fit a 2D polynomial
    in v & J to get the frequency offset
    """
    m_init = modeling.polynomial.Polynomial2D(2)
    fit = modeling.fitting.LevMarLSQFitter()
    xx, yy = species_diff_table['vu'], species_diff_table['Ju']
    zz = species_diff_table['Splat-Barton']
    model_fit = fit(m_init, xx, yy, zz)
    return model_fit

def match_splat_barton(splat, barton):
    closest_inds = []
    for row in splat:
        closest = np.argmin(np.abs(barton['Freq']-row['Freq']))
        closest_inds.append(closest)

    barton_ = barton[np.array(closest_inds)]
    result = table.hstack([splat, barton_])
    result.add_column(table.Column(data=result['Freq_1'] - result['Freq_2'], name='Splat-Barton'))
    return result

KCls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' KCl', line_lists=['SLAIM']))
KCl = load_barton('39K-35Cl')

KCl_diff = match_splat_barton(KCls, KCl[KCl['vu']<4])
KCl_offset_model = get_offset_model(KCl_diff)

KCl['Freq'] = KCl['Freq'] + KCl_offset_model(KCl['vu'], KCl['Ju'])



K37Cls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' K37Cl', line_lists=['SLAIM']))
K37Cl = load_barton('39K-37Cl')

K37Cl_diff = match_splat_barton(K37Cls, K37Cl[K37Cl['vu']<4])
K37Cl_offset_model = get_offset_model(K37Cl_diff)

K37Cl['Freq'] = K37Cl['Freq'] + K37Cl_offset_model(K37Cl['vu'], K37Cl['Ju'])


K41Cls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='41KCl', line_lists=['SLAIM']))
K41Cl = load_barton('41K-35Cl')

K41Cl_diff = match_splat_barton(K41Cls, K41Cl[K41Cl['vu']<4])
K41Cl_offset_model = get_offset_model(K41Cl_diff)

K41Cl['Freq'] = K41Cl['Freq'] + K41Cl_offset_model(K41Cl['vu'], K41Cl['Ju'])


K41Cl37 = load_barton('41K-37Cl')
# using 41KCl b/c there's nothing to fit here
K41Cl37['Freq'] = K41Cl37['Freq'] + K41Cl_offset_model(K41Cl37['vu'], K41Cl37['Ju'])

AlCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlCl'))
#AlF = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlF'))
#NaF = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaF'))
#NaO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaO'))
#NaOH = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaOH'))
#NaCH = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCH'))
NaCN = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCN'))
CaCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' CaCl'))
AlO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' AlO'))
#AlO = AlO[np.array([len(row['QNs']) < 10 for row in AlO])]
NaCls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCl', line_lists=['SLAIM']))
NaCl = load_barton('23Na-35Cl')

NaCl_diff = match_splat_barton(NaCls, NaCl[NaCl['vu']<4])
NaCl_offset_model = get_offset_model(NaCl_diff)

NaCl['Freq'] = NaCl['Freq'] + NaCl_offset_model(NaCl['vu'], NaCl['Ju'])

# use Cabezas instead of Barton
cabezasNaCl = Vizier(row_limit=1e9).get_catalogs('J/ApJ/825/150')[1]
cabezasNaCl['nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'] = u.Quantity(cabezasNaCl['nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'], u.cm**-1).to(u.MHz, u.spectral())
cabezasNaCl['x_nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'] = 'MHz'
cabezasNaCl.add_column(Column(name='Freq', data=u.Quantity(cabezasNaCl['nuCalc'], u.MHz).to(u.GHz)))
cabezasNaCl.rename_column('J1','Ju')
cabezasNaCl.rename_column('J0','Jl')
cabezasNaCl.rename_column('V1','vu')
cabezasNaCl.rename_column('V0','vl')
cabezasNaCl.rename_column('Eup','E_U')
NaCl = cabezasNaCl[cabezasNaCl['Iso'] == b'35']
NaCl.add_column(Column(name='Species',
                       data=['23Na-35Cl v={0}-{1} J={2}-{3}'
                             .format(row['vu'], row['vl'], row['Ju'], row['Jl'])
                             for row in NaCl]))


Na37Cls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' Na37Cl'))
Na37Cl = load_barton('23Na-37Cl')

Na37Cl_diff = match_splat_barton(Na37Cls, Na37Cl[Na37Cl['vu']<4])
Na37Cl_offset_model = get_offset_model(Na37Cl_diff)

Na37Cl['Freq'] = Na37Cl['Freq'] + Na37Cl_offset_model(Na37Cl['vu'], Na37Cl['Ju'])

Na37Cl = cabezasNaCl[cabezasNaCl['Iso'] == b'37']
Na37Cl.add_column(Column(name='Species',
                         data=['23Na-37Cl v={0}-{1} J={2}-{3}'
                               .format(row['vu'], row['vl'], row['Ju'], row['Jl'])
                               for row in NaCl]))


MgCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' MgCl'))
#MgCl = [row for row in MgCl if len(row['Resolved QNs']) < 20]
#not detected:
HCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' HCl'))
#NaCN = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCN')
#NaO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaO')
#FeCO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' FeCO')
#whiffs at 230 GHz:
# MgCCH = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='MgCCH')
#too many isotopologues don't have anything
# SiS = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='Silicon monosulfide')
#LOTS of lines.  No way.
# MnO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='MnO')
# SO is real, but only the main isotopologue?
#SO = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='Sulfur Monoxide')
SO = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' SO '))
S34O = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' 34SO '))
# only low-J lines of SO2...
# no, SO2 isn't really believable.
SO2 = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' SO2',
                                 energy_max=500, energy_type='eu_k'))

SiO = load_barton('28Si-16O')
SiO17 = load_barton('28Si-17O')
#Si29O = load_barton('29Si-16O')
Si30O = load_barton('30Si-16O')
sio_tables = [SiO,
              SiO17,
              #Si29O,
              Si30O]

SiS = load_barton('28Si-32S')
SiS33 = load_barton('28Si-33S')
Si30S = load_barton('30Si-32S')
Si29S = load_barton('29Si-32S')
sis_tables = [SiS,
              SiS33,
              Si29S,
              Si30S]


salt_colors = ['b', 'm', 'darkgreen', 'orange', 'c', 'y']
salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]

if __name__ == "__main__":

    # investigate difference between splatalogue & Barton cats
    KCl_diff = match_splat_barton(KCls, KCl[KCl['vu']<4])
    NaCl_diff = match_splat_barton(NaCls, NaCl[NaCl['vu']<6])

    import pylab as pl
    pl.figure(1).clf()
    pl.plot(KCl_diff['Freq_1'], KCl_diff['Splat-Barton'], '.')
    pl.title("KCl")
    pl.xlabel("Frequency (GHz)")
    pl.ylabel("Splatalogue - Barton frequency difference (GHz)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_KCl.pdf'))

    pl.clf()
    pl.plot(KCl_diff['Freq_1'], (KCl_diff['Splat-Barton'].quantity/KCl_diff['Freq_1'].quantity*constants.c).to(u.km/u.s), '.')
    pl.title("KCl")
    pl.xlabel("Frequency (GHz)")
    pl.ylabel("Splatalogue - Barton frequency difference (km/s)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_KCl_kms.pdf'))

    pl.clf()
    pl.plot(KCl_diff['EU_K'], (KCl_diff['Splat-Barton'].quantity/KCl_diff['Freq_1'].quantity*constants.c).to(u.km/u.s), '.')
    pl.title("KCl")
    pl.xlabel("E$_U$ [K]")
    pl.ylabel("Splatalogue - Barton frequency difference (km/s)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_KCl_kms_vs_EU.pdf'))


    pl.figure(2).clf()
    pl.plot(NaCl_diff['Freq_1'], NaCl_diff['Splat-Barton'], '.')
    pl.title("NaCl")
    pl.xlabel("Frequency (GHz)")
    pl.ylabel("Splatalogue - Barton frequency difference (GHz)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_NaCl.pdf'))

    pl.clf()
    pl.title("NaCl")
    pl.xlabel("Frequency (GHz)")
    pl.plot(NaCl_diff['Freq_1'], (NaCl_diff['Splat-Barton'].quantity/NaCl_diff['Freq_1'].quantity*constants.c).to(u.km/u.s), '.')
    pl.ylabel("Splatalogue - Barton frequency difference (km/s)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_NaCl_kms.pdf'))

    pl.clf()
    pl.title("NaCl")
    pl.xlabel("E$_U$ [K]")
    pl.plot(NaCl_diff['EU_K'], (NaCl_diff['Splat-Barton'].quantity/NaCl_diff['Freq_1'].quantity*constants.c).to(u.km/u.s), '.')
    pl.ylabel("Splatalogue - Barton frequency difference (km/s)")
    pl.savefig(paths.fpath('Splatalogue-Barton_comparison_NaCl_kms_vs_EU.pdf'))


    nacl_offset_model = get_offset_model(NaCl_diff)

    pl.subplot(2,1,1)
    pl.title("NaCl")
    pl.xlabel("vu")
    pl.plot(NaCl_diff['vu'], (NaCl_diff['Splat-Barton'].quantity), '.')
    pl.plot(NaCl_diff['vu'], nacl_offset_model(NaCl_diff['vu'], NaCl_diff['Ju']), '.')
    pl.subplot(2,1,2)
    pl.xlabel("Ju")
    pl.plot(NaCl_diff['Ju'], (NaCl_diff['Splat-Barton'].quantity), '.')
    pl.plot(NaCl_diff['Ju'], nacl_offset_model(NaCl_diff['vu'], NaCl_diff['Ju']), '.')
    pl.ylabel("Splatalogue - Barton frequency difference")

    kcl_offset_model = get_offset_model(KCl_diff)

    pl.figure(1).clf()
    pl.subplot(2,1,1)
    pl.title("KCl")
    pl.xlabel("vu")
    pl.plot(KCl_diff['vu'], (KCl_diff['Splat-Barton'].quantity), '.')
    pl.plot(KCl_diff['vu'], kcl_offset_model(KCl_diff['vu'], KCl_diff['Ju']), '.')
    pl.subplot(2,1,2)
    pl.xlabel("Ju")
    pl.plot(KCl_diff['Ju'], (KCl_diff['Splat-Barton'].quantity), '.')
    pl.plot(KCl_diff['Ju'], kcl_offset_model(KCl_diff['vu'], KCl_diff['Ju']), '.')
    pl.ylabel("Splatalogue - Barton frequency difference")

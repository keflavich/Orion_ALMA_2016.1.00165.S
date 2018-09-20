from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
import paths
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt

def load_barton(species):
    tbl = table.Table.read(paths.salty(species+"_rotational_transitions.ipac"),
                           format='ascii.ipac')
    tbl.add_column(table.Column(name='Species', data=[species+row['QNs'] for row in tbl]))
    tbl.rename_column('Frequency', 'Freq')
    tbl = tbl[(tbl['Freq']>10) & (tbl['Freq']<400) & (tbl['E_U'] < 1e4) & (tbl['vu']<=10)]
    return tbl

kcl_offset = (1+17*u.km/u.s/constants.c).decompose()
nacl_offset = (1+3.*u.km/u.s/constants.c).decompose()

KCls = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' KCl'))
KCl = load_barton('39K-35Cl')
KCl['Freq'] = KCl['Freq']*kcl_offset

#K37Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' K37Cl'))
K37Cl = load_barton('39K-37Cl')
K37Cl['Freq'] = K37Cl['Freq']*kcl_offset
#K41Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name='41KCl'))
K41Cl = load_barton('41K-35Cl')
K41Cl['Freq'] = K41Cl['Freq']*kcl_offset
K41Cl37 = load_barton('41K-37Cl')
K41Cl37['Freq'] = K41Cl37['Freq']*kcl_offset
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
#NaCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' NaCl'))
NaCl = load_barton('23Na-35Cl')
NaCl['Freq'] = NaCl['Freq']*nacl_offset
#Na37Cl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' Na37Cl'))
Na37Cl = load_barton('23Na-37Cl')
Na37Cl['Freq'] = Na37Cl['Freq']*nacl_offset

MgCl = mt(Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' MgCl'))
#MgCl = [row for row in MgCl if len(row['Resolved QNs']) < 20]
#not detected:
# HCl = Splatalogue.query_lines(80*u.GHz, 400*u.GHz, chemical_name=' HCl')
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


salt_colors = ['b', 'm', 'darkgreen', 'orange', 'c', 'y']
salt_tables = [KCl, K37Cl, K41Cl, NaCl, Na37Cl, K41Cl37]

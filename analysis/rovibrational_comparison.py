import numpy as np
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord
from astropy import modeling
from astropy.table import Column
import paths
from astroquery.splatalogue import Splatalogue
from astroquery.splatalogue.utils import minimize_table as mt
from astroquery.vizier import Vizier

def load_barton(species):
    tbl = table.Table.read(paths.salty(species+"_transitions.ipac"),
                           format='ascii.ipac')
    tbl.add_column(table.Column(name='Species', data=[species+row['QNs'] for row in tbl]))
    tbl.rename_column('Frequency', 'Freq')
    return tbl

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
    result.add_column(table.Column(data=(result['Freq_1'] - result['Freq_2'])/result['Freq_2']*3e5, name='Splat-Barton_kms'))
    return result

KCls = mt(Splatalogue.query_lines(8*u.GHz, 40000*u.GHz, chemical_name=' KCl', line_lists=['SLAIM']))
KCl = load_barton('39K-35Cl')

KCl_diff = match_splat_barton(KCls, KCl[KCl['vu']<4])
KCl_offset_model = get_offset_model(KCl_diff)


NaCls = mt(Splatalogue.query_lines(8*u.GHz, 40000*u.GHz, chemical_name=' NaCl', line_lists=['SLAIM']))
NaCl = load_barton('23Na-35Cl')

NaCl_diff = match_splat_barton(NaCls, NaCl[NaCl['vu']<4])
NaCl_offset_model = get_offset_model(NaCl_diff)

cabezasNaCl = Vizier(row_limit=1e9).get_catalogs('J/ApJ/825/150')[1]
cabezasNaCl['nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'] = u.Quantity(cabezasNaCl['nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'], u.cm**-1).to(u.MHz, u.spectral())
cabezasNaCl['x_nuCalc'][cabezasNaCl['x_nuCalc'] == 'cm-1'] = 'MHz'
cabezasNaCl.add_column(Column(name='Freq', data=u.Quantity(cabezasNaCl['nuCalc'], u.MHz).to(u.GHz)))

cabezas_match = match_splat_barton(cabezasNaCl, NaCl)
splat_cabezas_match = match_splat_barton(NaCls, cabezasNaCl)

import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.utils.console import ProgressBar
from astropy import log
import bz2
import paths

import requests


GHz_to_K = (1*u.GHz).to(u.eV, u.spectral()).to(u.K,
                                               u.temperature_energy()).value
cm_to_K = (1*u.cm**-1).to(u.eV, u.spectral()).to(u.K,
                                                 u.temperature_energy()).value

class ExoMol(object):
    def __init__(self, name, fullname, max_energy=15000*u.K, load_raw=False,
                 populate_transitions=True):

        self.base_url = 'http://www.exomol.com/{loc}/{name}/{fullname}/Barton/'.format(name=name,
                                                                                       fullname=fullname,
                                                                                       loc='{loc}')

        statesfile = paths.salty('{fullname}__Barton.states.bz2'.format(fullname=fullname))
        if os.path.exists(statesfile) and not load_raw:
            log.info("Loading state data from disk file {0}".format(statesfile))
            with open(statesfile, 'rb') as fh:
                content = fh.read()
            statesdata_unzip = bz2.BZ2Decompressor().decompress(content)
        else:
            log.info("Retrieving state data for {0}".format(fullname))
            statesdata = requests.get(self.base_url.format(loc='db')+'{fullname}__Barton.states.bz2'.format(fullname=fullname))
            statesdata_unzip = bz2.BZ2Decompressor().decompress(statesdata.content)
            with open(statesfile, 'wb') as fh:
                fh.write(statesdata.content)

        self.states = ascii.read(statesdata_unzip.decode().split('\n'))

        states_col_mapping = {'col1': ('state_id', None),
                              'col2': ('state_energy', u.cm**-1),
                              'col3': ('degeneracy', None),
                              'col4': ('J', None),
                              'col5': ('unknown', None),
                              'col6': ('v', None),
                             }
        for key,(name,unit) in states_col_mapping.items():
            self.states.rename_column(key, name)
            if unit is not None:
                self.states[name].unit = unit

        populated_transitions_fn = paths.salty('{fullname}_transitions.ipac'.format(fullname=fullname))
        rotpopulated_transitions_fn = paths.salty('{fullname}_rotational_transitions.ipac'.format(fullname=fullname))

        if os.path.exists(populated_transitions_fn) and not load_raw:
            log.info("Loading full transition table {0}".format(populated_transitions_fn))
            self.transitions = Table.read(populated_transitions_fn, format='ascii.ipac')
        else:

            transfile = paths.salty('{fullname}__Barton.trans.bz2'.format(fullname=fullname))
            if os.path.exists(transfile):
                log.info("Loading transition data from disk file {0}".format(transfile))
                with open(transfile, 'rb') as fh:
                    content = fh.read()
                transitionsdata_unzip = bz2.BZ2Decompressor().decompress(content)
            else:
                log.info("Retrieving transition data")
                transitionsdata = requests.get(self.base_url.format(loc='db')+'{fullname}__Barton.trans.bz2'.format(fullname=fullname))
                transitionsdata_unzip = bz2.BZ2Decompressor().decompress(transitionsdata.content)
                with open(transfile, 'wb') as fh:
                    fh.write(transitionsdata.content)

            self.transitions = ascii.read(transitionsdata_unzip.decode().split('\n'))

            transitions_col_mapping = {'col1': ('upper_state', None),
                                       'col2': ('lower_state', None),
                                       'col3': ('Aij', u.s**-1),
                                       'col4': ('Frequency_cm', u.cm**-1),
                                      }
            for key,(name,unit) in transitions_col_mapping.items():
                self.transitions.rename_column(key, name)
                if unit is not None:
                    self.transitions[name].unit = unit

            if max_energy is not None:
                max_energy = max_energy.to(u.K)
                state_energies = self.states['state_energy'] * cm_to_K

                selected_states = state_energies.data < max_energy.value
                # need to keep the indexing linear from zero
                max_state_id = len(selected_states) - np.argmax(selected_states[::-1])

                self.states = self.states[:max_state_id]

                sel_transitions = self.transitions['upper_state'] <= max_state_id
                log.info("Selecting down from {0} to {1} transitions"
                         .format(len(self.transitions), sel_transitions.sum()))
                self.transitions = self.transitions[sel_transitions]

            if populate_transitions:
                self.populate_table(max_energy)

            self.transitions.write(populated_transitions_fn,
                                   format='ascii.ipac',
                                  )

        if not os.path.exists(rotpopulated_transitions_fn):
            # rotational only
            rotational = self.transitions['vu'] == self.transitions['vl']
            self.transitions[rotational].write(rotpopulated_transitions_fn,
                                               format='ascii.ipac')

    def populate_table(self, max_energy):

        log.info("Populating Table w/E_U, E_L, frequency")
        self.transitions.add_column(Column(name='Frequency', unit=u.GHz,
                                           data=np.zeros(len(self.transitions)) + np.nan))
        self.transitions.add_column(Column(name='Wavelength', unit=u.um,
                                           data=np.zeros(len(self.transitions)) + np.nan))
        self.transitions.add_column(Column(name='E_U', unit=u.K,
                                           data=np.zeros(len(self.transitions)) + np.nan))
        self.transitions.add_column(Column(name='E_L', unit=u.K,
                                           data=np.zeros(len(self.transitions)) + np.nan))
        self.transitions.add_column(Column(name='QNs',
                                           data=np.empty(len(self.transitions), dtype='S25')))
        self.transitions.add_column(Column(name='vu',
                                           data=np.empty(len(self.transitions), dtype='int')))
        self.transitions.add_column(Column(name='vl',
                                           data=np.empty(len(self.transitions), dtype='int')))
        self.transitions.add_column(Column(name='Ju',
                                           data=np.empty(len(self.transitions), dtype='int')))
        self.transitions.add_column(Column(name='Jl',
                                           data=np.empty(len(self.transitions), dtype='int')))
        self.transitions.add_column(Column(name='gu',
                                           data=np.empty(len(self.transitions), dtype='int')))
        self.transitions.add_column(Column(name='gl',
                                           data=np.empty(len(self.transitions), dtype='int')))

        self.states.add_column(Column(name='state_frequency', unit=u.GHz,
                                      data=u.Quantity(self.states['state_energy'],
                                                      u.cm**-1).to(u.GHz,
                                                                   u.spectral())))

        max_energy = max_energy.to(u.K)

        for row in ProgressBar(self.transitions):
            lowerstate = row['lower_state']-1
            upperstate = row['upper_state']-1

            e_u = u.Quantity(self.states[upperstate]['state_frequency'], u.GHz).value
            row['E_U'] = e_u * GHz_to_K

            if e_u * GHz_to_K > max_energy.value:
                continue

            e_l = u.Quantity(self.states[lowerstate]['state_frequency'], u.GHz).value
            frequency = e_u - e_l

            row['Frequency'] = frequency
            row['E_L'] = e_l * GHz_to_K

            row['QNs'] = ('v={0}-{1} J={2}-{3}'
                          .format(int(self.states[upperstate]['v']),
                                  int(self.states[lowerstate]['v']),
                                  int(self.states[upperstate]['J']),
                                  int(self.states[lowerstate]['J']),))
            row['Ju'] = int(self.states[upperstate]['J'])
            row['Jl'] = int(self.states[lowerstate]['J'])
            row['vu'] = int(self.states[upperstate]['v'])
            row['vl'] = int(self.states[lowerstate]['v'])
            row['gu'] = int(self.states[upperstate]['degeneracy'])
            row['gl'] = int(self.states[lowerstate]['degeneracy'])

        self.transitions['Wavelength'] = self.transitions['Frequency'].quantity.to(u.um, u.spectral())
        self.transitions = self.transitions[self.transitions['E_U'] < max_energy.value]


if __name__ == "__main__":

    #k39cl35_exo = ExoMol('KCl', '39K-35Cl')
    #k39cl35_exo = ExoMol('KCl', '39K-35Cl', max_energy=None, load_raw=True,
    #                     populate_transitions=False)
    na23cl35_exo = ExoMol('NaCl', '23Na-35Cl', max_energy=100000*u.K, load_raw=True)
    #na23cl37_exo = ExoMol('NaCl', '23Na-37Cl')
    #k41cl35_exo = ExoMol('KCl', '41K-35Cl')
    #k41cl37_exo = ExoMol('KCl', '41K-37Cl')
    #k39cl37_exo = ExoMol('KCl', '39K-37Cl')

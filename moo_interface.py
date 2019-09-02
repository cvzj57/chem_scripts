import os
import sys
import json
import argparse
import multivariate_optimisation

empty_setup_file = {
    '_comment': 'Settings/control file for MOO.'
                '1) Specify ECPs to modify with the ecp_locators key '
                '2) Specify any pseudo geometry to modify in pseudo_geometry key '
                '3) Switch on different optimisation criteria with the optimise_* keys '
                '  and supply required reference data in tracked_* keys '
                '4) Specify starting variables yourself (supply_own_starting_guesses) or give a number '
                '  of semi-random (Gaussian-weighted) seeds to generate '
                '5) run the multivariate_optimisation.py script '
                '6) Good luck.',
    'calc_folder_path': '.',
    'optimise_pseudo_geometry': True,
    'pseudo_geometry_type': '',
    'optimise_with_orbitals': True,
    "optimise_with_total_gaps": False,
    'seeded_optimisation': False,
    'pseudo_geometry': {
        'indices_of_pseudo_carbons': [1],
        'potential_set_distance_guess': 0.5,
        'potential_set_split_distance_guess': 0.25,
        'potential_set_distance_bounds': (0.5, 1.0),
        'potential_set_split_distance_bounds': (0.25, 1.0)
    },
    'ecp_locators': [
        {'line_type': '$ecp',
         'basis_ecp_name': 'c ecp-p',
         'orbital_descriptor': 'p-f',
         },
        {'line_type': '$ecp',
         'basis_ecp_name': 'he ecp-s',
         'orbital_descriptor': 's-f',
         }
    ],
    'tracked_orbitals': [{
        'irrep': '1a',
        'spin': 'mos',
        'reference_energy': 0.0,
        'orbital_error_multiplier': 1.0
    },
    ],
    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)
    'initial_guesses': [0.1, 0.1, 0.1, 0.1],
    'initial_seed_number': 1,
    "tracked_total_gaps": [
        {
            "reference_gap_energy": 3.53295493,
            "gap_folder_path": "triplet"
        },
        {
            "reference_gap_energy": 9.09111851,
            "gap_folder_path": "cation"
        }
    ]
}


class MOOInterface:
    """Text interface for setting up MOO optimisations."""
    def __init__(self):
        # generic setup stuff
        self.setup_file_name = 'opt.moo'
        self.optdata = empty_setup_file
        self.line_separator = '------------------------'
        self.potential_reference = '*** [1] Punter, Alexander and Nava, Paola and Carissan, Yannick. Int. J. Quantum Chem. 2019. e25914.'
        self.menu_actions = {
            'moo_file_check1': 'confirm_opt_run',
            'moo_file_check2': 'initial_menu',
            'initial_menu': 'initial_menu',
            'initial_menu1': 'potentialisation_menu',
            'initial_menu2': 'optimisation_menu',
            'initial_menub': 'initial_menu',
            'initial_menuq': 'exit',
            # Optimisation menus
            'optimisation_menu': 'optimisation_menu',
            'optimisation_menu1': 'ecp_locator_menu',
            'optimisation_menu2': 'mo_criterion_menu',
            'optimisation_menu3': 'optimisation_menu',
            'optimisation_menu4': 'seed_options_menu',
            'optimisation_menu5': 'geomopt_options_menu',
            'optimisation_menu6': 'confirm_run_opt',
            'optimisation_menub': 'initial_menu',
            ## ECP locators
            'ecp_locator_menu1': 'ecp_locator_menu',
            'ecp_locator_menu2': 'optimisation_menu',
            'ecp_locator_menub': 'optimisation_menu',
            ## Orbital Criteria
            'mo_criterion_menu1': 'mo_criterion_menu',
            'mo_criterion_menu2': 'optimisation_menu',
            'mo_criterion_menub': 'optimisation_menu',
            ## Seed options
            'seed_options_menu': 'seed_options_menu',
            'seed_options_menu1': 'seed_options_menu',
            'seed_options_menu2': 'seed_number_menu',
            'seed_options_menub': 'optimisation_menu',
            ## Geomopt options
            'geomopt_options_menu1': 'geomopt_options_menu',
            'geomopt_options_menu2': 'geomopt_options_menu',
            'geomopt_options_menu3': 'geomopt_options_menu',
            'geomopt_options_menu4': 'geomopt_options_menu',
            # Potentialisation menus
            'potentialisation_menu1': 'potentialisation_menu',
            'potentialisation_menu2': 'potentialisation_menu',
            'potentialisation_menu3': 'potentialisation_menu',
            'potentialisation_menu4': 'potentialisation_menu',
        }
        print('''
                    #######################################
                     The Multi-Orbital Optimiser   (___)
                                                  <(o o)>______
                              (or MOO)              ../ ` # # \`;   
                                                      \ ,___, /
                        A Punter, CTOM group           ||   ||  
                     Université d'Aix-Marseille        ^^   ^^  
                    #######################################
            ''')

    def initial_menu(self):
        print('What would you like to do?')
        print('1: Potentialise a molecule.')
        print('2: Optimise new potentials.')
        print("q: Quit")
        choice = input(" >>  ")
        self.exec_menu(choice, 'initial_menu')

    def optimisation_menu(self):
        self.write_setup_file()
        print(self.line_separator)
        print('OPTIMISATION MENU')
        self.showoptsettings()
        print('1: Add ECP functions')
        print('2: Add MO criterion')
        print('3: Add total energy difference criterion')
        print('4: Semi-random seed options')
        print('5: Potential geometry optimisation options')
        print('6: Run optimisation now')
        print('b: Go back')
        choice = input(" >>  ")
        self.exec_menu(choice, 'optimisation_menu')

    def confirm_run_opt(self):
        if input('Run optimisation (y/n)?') == 'y':
            optimiser = multivariate_optimisation.Optimiser()
            optimiser.run()
        else:
            self.gotomenu('initial_menu')

    def ecp_locator_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: ECP LOCATORS')
        print('Please specify ECPs for MOO to alter')
        print('ECP name (e.g. "c ecp-1"):')
        ecp_name = input(' >>  ')
        print('angular momentum (e.g. "p-f"):')
        ecp_l = input(' >>  ')
        self.add_ecp_locator(ecp_name, ecp_l)
        print('1. Add another ECP function.')
        print('2. Done adding ECP functions.')
        choice = input(" >>  ")
        self.exec_menu(choice, 'ecp_locator_menu')

    def mo_criterion_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: MO CRITERION')
        print('MO name (e.g. "1a\""):')
        irrep = input(' >>  ')
        print('reference energy (eV):')
        ref_E = input(' >>  ')
        print('spin: (alpha, beta, closed)')
        spin = input(' >>  ')
        self.add_mo_criterion(irrep, ref_E, spin)
        print('1. Add another ECP function.')
        print('2. Done adding ECP functions.')
        choice = input(" >>  ")
        self.exec_menu(choice, 'mo_criterion_menu')

    def total_gap_criterion_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: TOTAL E GAP CRITERION')
        print('Comparison calc folder name (e.g. "1st_ie"):')
        other_calc_folder = input(' >>  ')
        print('reference total E gap (eV):')
        ref_E = input(' >>  ')
        self.add_total_gap_criterion(other_calc_folder, ref_E)
        print('1. Add another ECP function.')
        print('2. Done adding ECP functions.')
        choice = input(" >>  ")
        self.exec_menu(choice, 'total_gap_criterion_menu')

    def seed_options_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: SEED OPTIONS')
        print('1. Turn on/off seeded optimisation. (currently %s)' % self.trueon(self.optdata['seeded_optimisation']))
        print('2. Set seed number.')
        print('b: Go back')
        choice = input(" >>  ")
        if choice == '1':
            self.toggle_seeded()
        self.exec_menu(choice, 'seed_options_menu')

    def seed_number_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: SEED NUMBER')
        print('Enter seed number:')
        choice = input(" >>  ")
        self.set_seed_number(choice)
        self.gotomenu('seed_options_menu')

    def geomopt_options_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: POTENTIAL GEOMETRY')
        print('If using non-atom centered potentials in the style of [1], MOO can optimise their geometry.')
        print('1. Turn on/off geometry optimisation. (currently %s)'
              % self.trueon(self.optdata['optimise_pseudo_geometry']))
        print('2. Add pseudocarbon indices.')
        print('3. Set potential set distance boundaries (α,β,γ,δ potentials). Currently: %s'
              % self.optdata['pseudo_geometry']['potential_set_distance_bounds'])
        print('4. Set potential set split distance boundaries (α,β potentials). Currently: %s'
              % self.optdata['pseudo_geometry']['potential_set_split_distance_bounds'])
        print('5. Set pseudogeometry type (currently: %s).' % self.optdata['pseudo_geometry_type'])
        print(self.potential_reference)
        choice = input(" >>  ")
        if choice == '1':
            self.toggle_geomopt()
        self.exec_menu(choice, 'geomopt_options_menu')

    def moo_file_check(self):
        opt_data_path = os.path.join(os.getcwd(), self.setup_file_name)
        if not os.path.isfile(opt_data_path):
            print('No MOO file, running setup...')
            self.optdata = empty_setup_file
            getattr(self, self.menu_actions['initial_menu'])()
        else:
            print('MOO file found.')
            with open(opt_data_path, 'r') as opt_data_file:
                self.optdata = json.load(opt_data_file)
            opt_data_file.close()
            print('Run optimisation?')
            print('1: Start MOO optimisation run.')
            print('2: Enter MOO setup.')
            choice = input(" >>  ")
            self.exec_menu(choice, menu='moo_file_check')

    def toggle_seeded(self):
        self.optdata['seeded_optimisation'] = not self.optdata['seeded_optimisation']

    def toggle_geomopt(self):
        self.optdata['optimise_pseudo_geometry'] = not self.optdata['optimise_pseudo_geometry']

    def add_ecp_locator(self, ecp_name, ecp_l):
        self.optdata['ecp_locators'].append(
            {'line_type': '$ecp',
             'basis_ecp_name': ecp_name,
             'orbital_descriptor': ecp_l+'-f',
             'functions_list': [{'coefficient': 0,
                                 'r^n': 1,
                                 'exponent': 0}]})
        return

    def add_mo_criterion(self, irrep, ref_E, spin):
        if spin == 'closed':
            spin = 'mos'
        self.optdata['tracked_orbitals'].append(
            {'irrep': irrep,
             'spin': spin,
             'reference_energy': ref_E})
        return

    def add_total_gap_criterion(self, other_calc_folder, ref_E):
        self.optdata['tracked_total_gaps'].append({
            "reference_gap_energy": ref_E,
            "gap_folder_path": other_calc_folder})
        return

    def set_seed_number(self, initial_seed_number):
        self.optdata['initial_seed_number'] = int(initial_seed_number)

    def potentialisation_menu(self):
        print(self.line_separator)
        print('POTENTIAL PLACEMENT MENU')
        print('Use this menu to place potentials in the manner of [1].')
        print('Beware! I cannot undo mistakes! Save your geometries before attempting this!')
        print('1. Guess for indices (guesses potentialisation including specified carbon atoms).')
        print('2. Guess for indices except (guesses potentialisation excluding specified carbon atoms).')
        print('2. Place sp2 non-atom-centered potentials.')
        print('3. Place sp3 non-atom-centered potentials.')
        print(self.potential_reference)
        choice = input(" >>  ")
        self.exec_menu(choice, menu='potentialisation_menu')

    def incl_guess_menu(self):
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: INCLUSIVE MOO GUESS')
        print('Specify indices of carbons to guess in Turbomole format (e.g. 1-10,12,16)')
        include = input(' >> ')
        # call guess incl
        print('Guess complete. You will need to specify electron occupation manually in the control file.')

    def excl_guess_menu(self):
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: EXCLUSIVE MOO GUESS')
        print('Specify indices of carbons to ignore in Turbomole format (e.g. 1-10,12,16). I\'ll guess the rest.')
        exclude = input(' >> ')
        # call guess excl
        print('Guess complete. You will need to specify electron occupation manually in the control file.')

    def pot_sp2_menu(self):
        pass

    def pot_sp3_menu(self):
        pass

    def trueon(self, boolean):
        if boolean is True:
            message = 'ON'
        else:
            message = 'OFF'
        return message

    def exec_menu(self, choice, menu=''):
        ch = choice.lower()
        if ch == '':
            getattr(self, self.menu_actions[menu])()
        else:
            action_id = menu + choice
            print(action_id)
            try:
                getattr(self, self.menu_actions[action_id])()
            except KeyError:
                print("Invalid selection, please try again.\n")
                getattr(self, self.menu_actions['initial_menu'])()
        return

    def gotomenu(self, menu):
        getattr(self, self.menu_actions[menu])()

    def showoptsettings(self):
        print(self.line_separator)
        print('CURRENT SETTINGS:')
        print('ECPs to optimise: ', [('%s' % option['basis_ecp_name'] + ' (%s)'
                                      % option['orbital_descriptor']) for option in self.optdata['ecp_locators']])
        print('Orbital Optimisation: ', self.trueon(self.optdata['optimise_with_orbitals']), '(orbitals: %s)'
              % [option['irrep'] for option in self.optdata['tracked_orbitals']])
        print('Total Gap Optimisation: ', self.trueon(self.optdata['optimise_with_total_gaps']), '(comparison folders: %s)'
              % [option['gap_folder_path'] for option in self.optdata['tracked_total_gaps']])
        print('Seeded Optimisation: ', self.trueon(self.optdata['seeded_optimisation']), '(seeds: %s)'
              % self.optdata['initial_seed_number'])
        print('Geometry Optimisation: ', self.trueon(self.optdata['optimise_pseudo_geometry']), '(%s carbons: %s)'
              % (self.optdata['pseudo_geometry_type'], self.optdata['pseudo_geometry']['indices_of_pseudo_carbons']))
        print(self.line_separator)

    # Back to main menu
    def back(self):
        getattr(self, self.menu_actions['main_menu'])()

    # Exit program
    def exit(self):
        sys.exit()

    def write_setup_file(self):
        with open(os.path.join(self.optdata['calc_folder_path'], self.setup_file_name), 'w') as outfile:
            json.dump(self.optdata, outfile, sort_keys=True, indent=2)
        outfile.close()


if __name__ == "__main__":
    interface = MOOInterface()
    interface.moo_file_check()
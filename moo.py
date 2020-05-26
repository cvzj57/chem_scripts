import os
import sys
import json
import argparse
import multivariate_optimisation
import chem_lib.coord as coord

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
    'optimise_pseudo_geometry': False,
    'pseudo_geometry_type': 'sp2',
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
    ],
    'tracked_orbitals': [
    ],
    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)
    'initial_guesses': [],
    'initial_seed_number': 1,
    "tracked_total_gaps": [
    ]
}


class MOOInterface:
    """Text interface for setting up MOO optimisations."""
    def __init__(self):
        # generic setup stuff
        self.setup_file_name = 'opt.moo'
        self.optdata = empty_setup_file
        self.line_separator = '--------------------------------------------'
        self.references = '''[1] Punter, Alexander and Nava, Paola and Carissan, Yannick. Int. J. Quantum Chem. 2019. e25914.'''
        self.moo_logo = '''
                    #######################################
                     The Multi-Orbital Optimiser   (___)
                                                  <(o o)>______
                              (or MOO)              ../ ` # # \`;   
                                                      \ ,___, /
                        A Punter, CTOM group           ||   ||  
                      Aix-Marseille Université         ^^   ^^  
                    #######################################
            '''
        print(self.moo_logo)
        self.menu_actions = {
            'moo_file_check1': 'confirm_run_opt',
            'moo_file_check2': 'initial_menu',
            'initial_menu': 'initial_menu',
            'initial_menu1': 'potentialisation_menu',
            'initial_menu2': 'optimisation_menu',
            'initial_menu3': 'initial_menu',
            'initial_menub': 'initial_menu',
            'initial_menuq': 'exit',
            # Optimisation menus
            'optimisation_menu': 'optimisation_menu',
            'optimisation_menu1': 'ecp_locator_menu',
            'optimisation_menu2': 'optimisation_menu',
            'optimisation_menu3': 'mo_criterion_menu',
            'optimisation_menu4': 'optimisation_menu',
            'optimisation_menu5': 'total_gap_criterion_menu',
            'optimisation_menu6': 'seed_options_menu',
            'optimisation_menu7': 'geomopt_options_menu',
            'optimisation_menu8': 'initial_guess_menu',
            'optimisation_menu9': 'confirm_run_opt',
            'optimisation_menuh': 'optimisation_menu',
            'optimisation_menub': 'initial_menu',
            ## ECP locators
            'ecp_locator_menu': 'ecp_locator_menu',
            'ecp_locator_menu1': 'ecp_locator_menu',
            'ecp_locator_menu2': 'optimisation_menu',
            'ecp_locator_menub': 'optimisation_menu',
            ## Orbital Criteria
            'mo_criterion_menu': 'mo_criterion_menu',
            'mo_criterion_menu1': 'mo_criterion_menu',
            'mo_criterion_menu2': 'optimisation_menu',
            'mo_criterion_menub': 'optimisation_menu',
            ## Seed options
            'seed_options_menu': 'seed_options_menu',
            'seed_options_menu1': 'seed_options_menu',
            'seed_options_menu2': 'seed_number_menu',
            'seed_options_menub': 'optimisation_menu',
            ## Geomopt options
            'geomopt_options_menu': 'geomopt_options_menu',
            'geomopt_options_menu1': 'geomopt_options_menu',
            'geomopt_options_menu2': 'geom_indices_menu',
            'geomopt_options_menu3': 'geomopt_options_menu',
            'geomopt_options_menu4': 'geomopt_options_menu',
            'geomopt_options_menu5': 'geomopt_options_menu',
            'geomopt_options_menub': 'optimisation_menu',
            # Potentialisation menus
            'potentialisation_menu': 'potentialisation_menu',
            'potentialisation_menu1': 'potentialisation_menu',
            'potentialisation_menu2': 'incl_guess_menu',
            'potentialisation_menu3': 'excl_guess_menu',
            'potentialisation_menu4': 'pot_sp2_menu',
            'potentialisation_menu5': 'pot_sp3_menu',
            'potentialisation_menu6': 'repot_sp2_menu',
            'potentialisation_menu7': 'repot_sp3_menu',
            'potentialisation_menub': 'initial_menu',
        }
        self.coordcontrol = coord.CoordControl()
        self.post_guess_message = '''
        Guess complete. You will need to specify electron occupation manually 
        in the control file. Check the pseudification.log to see which atoms I replaced.'''
        self.optimisation_guide = '''
        OPTIMISATION GUIDE
        
        In order to optimise potentials with MOO, you will need control and coord files already set up with 
        pseudopotentials. 
        
        After this, you'll need to supply the ecp name and angular momentum of the potential to MOO using 
        the 'Add ECP functions' option. This lets MOO find the ecp to optimise in the basis file.
        
        Next, you will need either to supply a starting coefficient and exponent for each ECP using the 'Set initial 
        guesses...' option, or use the semi-random seed options. When supplying initial guesses, bear in mind that 
        smaller numbers tend to work better, and that exponents cannot be negative.
        
        Using the semi-random seed option means that MOO will run the optimisation a specified number of times using 
        semi-random guesses for the starting potential criteria. Several hundred seeds may be needed to be sure of 
        finding the best result within the boundary conditions, so this is much slower than a good starting guess.
        
        Now you should be ready to add optimisation criteria. MOO currently has two types of criteria to use. These are:
        (1) Molecular Orbital energies, which are supplied with the 'Molecular Orbital energy criteria').
        (2) Total energy differences, in which MOO will try to fit the difference between two different pseudopotential 
        calculations, both of which will need to be prepared in advance. Use the 'Total energy difference criteria'
        options to add these.
        
        Finally, the 'Potential geometry optimisation options' are to be used only with pseudofragments as described in
        Reference [1]. With this option the positions of non-atom-centered pseudopotentials will be adjusted by MOO
        during the optimisation. You will need to specify the indices of pseudocarbon atoms and the type of 
        hybridisation.
        '''
        try:
            self.coordcontrol.read_coords()
        except Exception as e:
            print('I can\'t see a coord file in this directory.')
            self.exit()

    def initial_menu(self):
        print('What would you like to do?')
        print('1: Potentialise a molecule.')
        print('2: Optimise new potentials.')
        print('3. View references.')
        print("q: Quit")
        choice = input(" >>  ")
        if choice == '3':
            self.print_references()
        self.exec_menu(choice, 'initial_menu')

    def optimisation_menu(self):
        self.write_setup_file()
        print(self.line_separator)
        print('OPTIMISATION MENU')
        self.showoptsettings()
        print('1. Add ECP functions')
        print('2. Toggle MO energy optimisation (currently %s)'
              % self.trueon(self.optdata['optimise_with_orbitals']))
        print('3. Molecular Orbital energy criteria')
        print('4. Toggle total energy difference optimisation (currently %s)'
              % self.trueon(self.optdata['optimise_with_total_gaps']))
        print('5. Total energy difference criteria')
        print('6. Semi-random seed options')
        print('7. Potential geometry optimisation options')
        print('8. Set initial guesses for ECP parameters (only needed if not using semi-random seeds).')
        print('9. Run optimisation now')
        print('h. Help. Show optimisation guide.')
        print('b: Go back')
        choice = input(" >>  ")
        if choice == '2':
            self.toggle_mo_crit()
        elif choice == '4':
            self.toggle_totalgap_crit()
        elif choice == 'h':
            self.print_optimisation_guide()
        self.exec_menu(choice, 'optimisation_menu')

    def print_optimisation_guide(self):
        print(self.line_separator)
        print(self.optimisation_guide)

    def confirm_run_opt(self):
        print('Run optimisation (y/n)?')
        choice = input(' >> ')
        if choice == 'y':
            try:
                optimiser = multivariate_optimisation.Optimiser()
                optimiser.run()
            except Exception as e:
                print('Error: %s' % e)
        else:
            self.gotomenu('initial_menu')

    def ecp_locator_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: ECP LOCATORS')
        print('Please specify ECPs for MOO to optimise')
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
        ref_E = float(input(' >>  '))
        print('spin: (alpha, beta, closed)')
        spin = input(' >>  ')
        self.add_mo_criterion(irrep, ref_E, spin)
        print('1. Add another MO criterion.')
        print('2. Done adding MO criteria.')
        choice = input(" >>  ")
        self.exec_menu(choice, 'mo_criterion_menu')

    def total_gap_criterion_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: TOTAL E GAP CRITERION')
        print('Comparison calc folder name (e.g. "1st_ie"):')
        other_calc_folder = input(' >>  ')
        print('reference total E gap (eV):')
        ref_E = float(input(' >>  '))
        self.add_total_gap_criterion(other_calc_folder, ref_E)
        print('1. Add another Total E Gap criterion.')
        print('2. Done adding Total E Gap criteria.')
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
        print('If using non-atom centered potentials in the style of reference [1], MOO can optimise their geometry.')
        print('1. Turn on/off geometry optimisation. (currently %s)'
              % self.trueon(self.optdata['optimise_pseudo_geometry']))
        print('2. Add pseudocarbon indices.')
        print('3. Toggle pseudogeometry type (currently: %s).' % self.optdata['pseudo_geometry_type'])
        print('b: Go back')
        #print('3. Set potential set distance boundaries (α,β,γ,δ potentials). Currently: %s'
        #      % self.optdata['pseudo_geometry']['potential_set_distance_bounds'])
        #print('4. Set potential set split distance boundaries (α,β potentials). Currently: %s'
        #      % self.optdata['pseudo_geometry']['potential_set_split_distance_bounds'])
        choice = input(" >>  ")
        if choice == '1':
            self.toggle_geomopt()
        elif choice == '3':
            self.toggle_geomtype()
        self.exec_menu(choice, 'geomopt_options_menu')

    def geom_indices_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: PSEUDOCARBON INDICES')
        print('Enter indices of pseudocarbons for potential geometry optimisation in Turbomole format (e.g. 1-10,12,16)')
        choice = input(" >>  ")
        self.set_geom_indices(choice)
        self.gotomenu('geomopt_options_menu')

    def initial_guess_menu(self):
        print(self.line_separator)
        print('OPTIMISATION: INITIAL GUESSES')
        print('If not using semi-random seeds, please supply an initial coefficient and exponent for each ECP function.')
        initial_guesses = []
        for ecp_locator in self.optdata['ecp_locators']:
            print('Initial coefficient for %s, %s:'
                  % (ecp_locator['basis_ecp_name'], ecp_locator['orbital_descriptor']))
            coefficient = input()
            initial_guesses.append(float(coefficient))
            print('Initial exponent for %s, %s:'
                  % (ecp_locator['basis_ecp_name'], ecp_locator['orbital_descriptor']))
            exponent = input()
            initial_guesses.append(float(exponent))
        self.optdata['initial_guesses'] = initial_guesses
        self.gotomenu('optimisation_menu')

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

    def toggle_mo_crit(self):
        self.optdata['optimise_with_orbitals'] = not self.optdata['optimise_with_orbitals']

    def toggle_totalgap_crit(self):
        self.optdata['optimise_with_total_gaps'] = not self.optdata['optimise_with_total_gaps']

    def toggle_seeded(self):
        self.optdata['seeded_optimisation'] = not self.optdata['seeded_optimisation']

    def toggle_geomopt(self):
        self.optdata['optimise_pseudo_geometry'] = not self.optdata['optimise_pseudo_geometry']

    def toggle_geomtype(self):
        print(self.optdata['pseudo_geometry_type'])
        if self.optdata['pseudo_geometry_type'] == 'sp2':
            self.optdata['pseudo_geometry_type'] = 'sp3'
        elif self.optdata['pseudo_geometry_type'] == 'sp3':
            self.optdata['pseudo_geometry_type'] = 'sp2'

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
             'reference_energy': ref_E,
             'orbital_error_multiplier': 1.0})
        return

    def add_total_gap_criterion(self, other_calc_folder, ref_E):
        self.optdata['tracked_total_gaps'].append({
            "reference_gap_energy": ref_E,
            "gap_folder_path": other_calc_folder})
        return

    def set_seed_number(self, initial_seed_number):
        self.optdata['initial_seed_number'] = int(initial_seed_number)

    def set_geom_indices(self, indices_string):
        self.optdata['pseudo_geometry']['indices_of_pseudo_carbons'] = self.parse_coord_list(indices_string)

    def potentialisation_menu(self):
        print(self.line_separator)
        print('POTENTIAL PLACEMENT MENU')
        print('Use this menu to place potentials in the manner of reference [1].')
        print('WARNING: Beware! I cannot undo mistakes! Save your geometries before attempting this!')
        print('1. Guess potentialisation of whole molecule.')
        print('2. Guess for indices (guesses potentialisation including specified carbon atoms).')
        print('3. Guess for indices except (guesses potentialisation excluding specified carbon atoms).')
        print('4. Place sp2 non-atom-centered potentials.')
        print('5. Place sp3 non-atom-centered potentials.')
        print('6. Reposition sp2 non-atom-centered potentials.')
        print('7. Reposition sp3 non-atom-centered potentials.')
        print('b: Go back')
        choice = input(" >>  ")
        if choice == '1':
            self.reset_coords()
            self.coordcontrol.guess_potentialisation('')
            print(self.post_guess_message)
            self.gotomenu('initial_menu')
        else:
            self.exec_menu(choice, menu='potentialisation_menu')

    def incl_guess_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: INCLUSIVE MOO GUESS')
        print('Specify indices of carbons to guess in Turbomole format (e.g. 1-10,12,16)')
        include = input(' >> ')
        self.coordcontrol.guess_potentialisation(['guess', 'incl', include])
        print(self.post_guess_message)
        self.gotomenu('initial_menu')

    def excl_guess_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: EXCLUSIVE MOO GUESS')
        print('Specify indices of carbons to ignore in Turbomole format (e.g. 1-10,12,16). I\'ll guess the rest.')
        exclude = input(' >> ')
        self.coordcontrol.guess_potentialisation(['guess', 'excl', exclude])
        print(self.post_guess_message)
        self.gotomenu('initial_menu')

    def pot_sp2_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: POTENTIALISE SP2 CARBONS')
        print('Specify pseudocarbon indices in Turbomole format (e.g. 1-10,12,16). Leave blank for all.')
        indices_string = input(" >>  ")
        # indices = self.parse_coord_list(indices_string)
        self.coordcontrol.pseudopotentialise_molecule('sp2', indices_string)
        print('Atoms pseudopotentialised.')
        self.gotomenu('potentialisation_menu')

    def pot_sp3_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: POTENTIALISE SP3 CARBONS')
        print('Specify sp3 atoms (e.g. hydrogen for methyl) indices in Turbomole format (e.g. 1-10,12,16).')
        indices_string = input(" >>  ")
        self.coordcontrol.pseudopotentialise_ethane_like_molecule(['', 'sp3', indices_string])
        print('Atoms pseudopotentialised.')
        self.gotomenu('potentialisation_menu')

    def repot_sp2_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: REPOTENTIALISE SP2 CARBONS')
        print('Specify pseudocarbon indices in Turbomole format (e.g. 1-10,12,16)')
        indices_string = input(" >>  ")
        indices = self.parse_coord_list(indices_string)
        print('Specify set distance in atomic units a.u.')
        set_dist = input(' >> ')
        print('Specify split distance in atomic units a.u.')
        split_dist = input(' >> ')
        for index in indices:
            self.coordcontrol.repseudopotentialise_sp2_atom(index, set_dist, split_dist)
        print('New pseudo distances set.')
        self.gotomenu('potentialisation_menu')

    def repot_sp3_menu(self):
        self.reset_coords()
        print(self.line_separator)
        print('POTENTIAL PLACEMENT: REPOTENTIALISE SP3 CARBONS')
        print('Specify pseudocarbon indices in Turbomole format (e.g. 1-10,12,16)')
        indices_string = input(" >>  ")
        indices = self.parse_coord_list(indices_string)
        print('Specify set distance in atomic units a.u.')
        set_dist = input(' >> ')
        for index in indices:
            self.coordcontrol.set_potential_distance_to(index, set_dist)
        print('New pseudo distances set.')
        self.gotomenu('potentialisation_menu')

    def trueon(self, boolean):
        if boolean is True:
            message = 'ON'
        else:
            message = 'OFF'
        return message

    def exec_menu(self, choice, menu=''):
        action_id = menu + choice
        try:
            getattr(self, self.menu_actions[action_id])()
        except KeyError:
            print("Invalid selection.\n")
            try:
                self.gotomenu(menu)
            except Exception as e:
                print(e)
                getattr(self, self.menu_actions['initial_menu'])()
        return

    def parse_coord_list(self, raw_input):
        print("Raw input received: %s" % raw_input)
        splitted = raw_input.split(',')
        final_coord_list = []
        for split_input in splitted:
            if '-' in split_input:
                input_range = split_input.split('-')
                final_coord_list.extend([*range(int(input_range[0]), int(input_range[1])+1)])
            else:
                final_coord_list.append(int(split_input))
        return final_coord_list

    def gotomenu(self, menu):
        getattr(self, self.menu_actions[menu])()

    def showoptsettings(self):
        print(self.line_separator)
        print('CURRENT SETTINGS:')
        print('ECPs to optimise:       ', [('%s' % option['basis_ecp_name'] + ' (%s)'
                                      % option['orbital_descriptor']) for option in self.optdata['ecp_locators']])
        print('Orbital Optimisation:   ', self.trueon(self.optdata['optimise_with_orbitals']), '(orbitals: %s)'
              % [option['irrep'] for option in self.optdata['tracked_orbitals']])
        print('Total Gap Optimisation: ', self.trueon(self.optdata['optimise_with_total_gaps']), '(comparison folders: %s)'
              % [option['gap_folder_path'] for option in self.optdata['tracked_total_gaps']])
        print('Seeded Optimisation:    ', self.trueon(self.optdata['seeded_optimisation']), '(seeds: %s)'
              % self.optdata['initial_seed_number'])
        print('Geometry Optimisation:  ', self.trueon(self.optdata['optimise_pseudo_geometry']), '(carbons of type %s, indices: %s)'
              % (self.optdata['pseudo_geometry_type'], self.optdata['pseudo_geometry']['indices_of_pseudo_carbons']))
        print('Initial guesses:        ', self.optdata['initial_guesses'])
        if not self.optdata['seeded_optimisation']:
            print('WARNING: As seeded optimisation is OFF, you will need to supply your own initial guess '
                  'for potential parameters.')
        print(self.line_separator)

    def print_references(self):
        print(self.line_separator)
        print('REFERENCES:')
        print(self.references)
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

    def reset_coords(self):
        self.coordcontrol.coord_list = []
        self.coordcontrol.read_coords()

if __name__ == "__main__":
    interface = MOOInterface()
    interface.moo_file_check()
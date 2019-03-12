#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.coord import CoordControl
from chem_lib.optimise import BasisControl
from chem_lib.gradient import GradientControl
from chem_lib.escf import ESCFControl
import scipy.optimize
import numpy
import json
import subprocess
import sys
import logging
import matplotlib.pyplot as plt

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
    'pseudo_geometry_type': '',
    'optimise_with_orbitals': True,
    'optimise_with_total_energy': False,
    'optimise_with_fixed_excitations': False,
    'optimise_with_flexible_excitations': False,
    'optimise_with_excitation_spectra': False,
    'optimise_with_homo_lumo_gap': False,
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
    'tracked_excitations': [{
        'repr': '1a',
        'reference_energy': 0.0,
        'reference_oscillation': 0.0,
        'reference_rotation': 0.0,
        'occ_orbitals': ['1a', '2a'],
        'reference_wavelength': 162.0,
        'electric_dipole_norm': 0.1
    },
    ],
    'tracked_homo_lumo_gap': 0.0,
    'tracked_total': 0.0,
    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)
    'initial_guesses': [0.1, 0.1, 0.1, 0.1],
    'seeded_optimisation': False,
    'initial_seed_number': 1
}

# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=os.path.join(os.getcwd(), 'moo.log'),
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)
# Now, we can log to the root logger, or any other logger. First the root...
logging.info('Logging initialised...')


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energies."""
    def __init__(self):
        # generic setup stuff
        self.setup_file_name = 'opt.moo'
        self.coord = CoordControl()
        self.basis = BasisControl()
        self.gradient = GradientControl()
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.calculation_type = 'dscf'
        self.calc_folder_path = '.'
        self.iterations = 0
        self.current_seed = 0
        self.moo_data_directory = 'moo_data'
        # optimisation specs
        self.ecp_locators = []

        self.optimise_with_total_energy = False
        self.tracked_total_energy = 0.0

        self.optimise_pseudo_geometry = False
        self.pseudo_geometry_type = ''
        self.pseudo_geometry = [] # geometric specifications for pseudo geometry optimisation

        self.optimise_with_orbitals = False
        self.tracked_orbitals = []  # locators for orbitals to attempt to optimise
        self.orbital_error_multiplier = 1.0

        self.optimise_with_fixed_excitations = False
        self.optimise_with_flexible_excitations = False
        self.optimise_with_excitation_spectra = False
        self.tracked_excitations = []  # locators for tddft excitations to attempt to optimise
        self.flexible_excitation_error_multiplier = 1.0
        self.flexible_excitation_osc_error_multiplier = 10.0
        self.excitation_spectra_error_multiplier = 1.0

        self.tracked_atom_gradients = []  # hashes of atoms from which to read gradient
        self.gradient_error_multiplier = 10.0  # multiplier gradient minimisation

        self.optimise_with_homo_lumo_gap = False
        self.tracked_homo_lumo_gap = 0.0
        self.homo_lumo_gap_error_multiplier = 1.0

        self.optimise_with_total_gaps = False
        self.tracked_total_gaps = []
        self.total_gaps_error_multiplier = 1.0

        self.seeded_optimisation = False
        self.initial_seed_number = 5
        self.initial_guesses = None
        with open(os.path.dirname(os.path.realpath(__file__))+'/chem_lib/optimiser_orbital_library.json', 'r') as json_file:
            self.orbital_library = json.load(json_file)


    def run(self):
        logging.info(
            '''
            #######################################
             The Multi-Orbital Optimiser   (___)
                                          '(o o)`______
                      (or MOO)              ../ ` # # \`~   
                                              \ ,___, /
                A Punter, CTOM group          //    //  
             Université d'Aix-Marseille       ^^    ^^  
            #######################################
            ''')
        opt_data_path = os.path.join(os.getcwd(), self.setup_file_name)
        if not os.path.isfile(opt_data_path):
            logging.info('No moo file, running setup...')
            self.setup_menu()
        else:
            logging.info('moo file found.')
        optdata = None
        try:
            with open(opt_data_path, 'r') as opt_data_file:
                optdata = json.load(opt_data_file)
                opt_data_file.close()
        except Exception as e:
            logging.info('I can\'t read your moo file. Check your JSON.\n %s' % e)
        self.optimise_pseudo_geometry = optdata.get('optimise_pseudo_geometry',
                                                    self.optimise_pseudo_geometry)
        self.pseudo_geometry_type = optdata.get('pseudo_geometry_type')
        self.optimise_with_orbitals = optdata.get('optimise_with_orbitals',
                                                  self.optimise_with_orbitals)
        self.optimise_with_total_energy = optdata.get('optimise_with_total_energy',
                                                      self.optimise_with_total_energy)
        self.optimise_with_fixed_excitations = optdata.get('optimise_with_fixed_excitations',
                                                           self.optimise_with_fixed_excitations)
        self.optimise_with_flexible_excitations = optdata.get('optimise_with_flexible_excitations',
                                                              self.optimise_with_flexible_excitations)
        self.optimise_with_excitation_spectra = optdata.get('optimise_with_excitation_spectra',
                                                            self.optimise_with_excitation_spectra)
        self.optimise_with_homo_lumo_gap = optdata.get('optimise_with_homo_lumo_gap',
                                                       self.optimise_with_homo_lumo_gap)
        self.optimise_with_total_gaps = optdata.get('optimise_with_total_gaps',
                                                   self.optimise_with_total_gaps)
        self.tracked_homo_lumo_gap = optdata.get('tracked_homo_lumo_gap')
        self.tracked_total_gaps = optdata.get('tracked_total_gaps', self.tracked_total_gaps)
        self.tracked_orbitals = optdata.get('tracked_orbitals')
        self.tracked_excitations = optdata.get('tracked_excitations')
        self.tracked_total_energy = optdata.get('tracked_total')
        self.calc_folder_path = optdata.get('calc_folder_path')
        self.ecp_locators = optdata.get('ecp_locators')
        self.pseudo_geometry = optdata.get('pseudo_geometry', {})
        self.initial_guesses = numpy.array(optdata['initial_guesses'])
        self.initial_seed_number = optdata.get('initial_seed_number', self.initial_seed_number)
        self.seeded_optimisation = optdata.get('seeded_optimisation', self.seeded_optimisation)
        if '-nosetup' in sys.argv:
            with open('errors.log', 'w') as error_file:
                try:
                    error_file.write('Starting run...')
                    self.run_multivariate_orbital_optimisation(numpy.array(self.initial_guesses))
                except Exception as e:
                    error_file.write('%s' % e)
            error_file.close()
        else:
            if input('Run now? y/n: ') == 'y':
                self.run_multivariate_orbital_optimisation(numpy.array(self.initial_guesses))

    def setup_menu(self):
        setup_file = empty_setup_file
        setup_file['calc_folder_path'] = input('Enter calc folder path: ')

        print('1. Enter tracked orbital key')
        print('2. Show orbital key library')
        print('3. Generate empty tracked orbital file')
        orbital_selection_made = False
        initial_guess_selection_made = False
        while orbital_selection_made is False:
            try:
                orbital_menu_choice = int(input('Enter choice: '))
                if orbital_menu_choice == 1:
                    orbital_key = input('Enter key: ')
                    try:
                        setup_file['tracked_orbitals'] = self.orbital_library['orbital_list'][orbital_key]
                        orbital_selection_made = True
                    except KeyError:
                        print('Not a real key')
                elif orbital_menu_choice == 2:
                    sorted_keys = sorted(list(self.orbital_library['orbital_list'].keys()))
                    for id_no, key in enumerate(sorted_keys):
                        print("%s: %s" % (id_no, sorted_keys[id_no]))
                    try:
                        chosen_id = int(input('Enter ID: '))
                        setup_file['tracked_orbitals'] = self.orbital_library['orbital_list'][sorted_keys[chosen_id]]
                    except Exception as error:
                        print('Not a real key... ', error)
                    orbital_selection_made = True
                elif orbital_menu_choice == 3:
                    orbital_selection_made = True
            except Exception as error:
                print("That's not a choice!", error)

        print('1. Enter initial guess key / #')
        print('2. Print initial guess library')
        print("3. Don't care (all values 0.1)")
        while initial_guess_selection_made is False:
            try:
                guess_menu_choice = int(input('Enter choice: '))
                if guess_menu_choice == 1:
                    guess_key = input('Enter key: ')
                    try:
                        setup_file['initial_guesses'] = self.orbital_library['initial_guess_list'][int(guess_key)]['guesses']
                    except KeyError:
                        print('Not a real key')
                    initial_guess_selection_made = True
                elif guess_menu_choice == 2:
                    for key, guess in enumerate(self.orbital_library['initial_guess_list']):
                        print("%s: %s \n    %s" % (key, guess['description'], guess['guesses']))
                    guess_key = input('Enter key: ')
                    try:
                        setup_file['initial_guesses'] = self.orbital_library['initial_guess_list'][int(guess_key)][
                            'guesses']
                    except KeyError:
                        print('Not a real key')
                    initial_guess_selection_made = True
                elif guess_menu_choice == 3:
                    initial_guess_selection_made = True
            except ValueError:
                print("That's not a choice!")
        if input('Add carbon s? y/n: ') == 'y':
            setup_file['ecp_locators'].append({'line_type': '$ecp',
                                               'basis_ecp_name': 'c ecp-p',
                                               'orbital_descriptor': 's-f',
                                               'functions_list': [{'coefficient': 0,
                                                                   'r^n': 1,
                                                                   'exponent': 0}]})

        with open(os.path.join(setup_file['calc_folder_path'], self.setup_file_name), 'w') as outfile:
            json.dump(setup_file, outfile, sort_keys=True, indent=2)
        outfile.close()
        print("opt.moo settings file created. Use this file to add optimisation "
              "of excitations, HOM-LUMO gap, total energy or gradients.")

    def create_empty_settings_file(self):
        setup_file = empty_setup_file
        with open(os.path.join(setup_file['calc_folder_path'], self.setup_file_name), 'w') as outfile:
            json.dump(setup_file, outfile, sort_keys=True, indent=2)
        outfile.close()
        print("opt.moo settings file created. Use this file to add optimisation "
              "of excitations, HOM-LUMO gap, total energy or gradients.")

    def remove_files(self, filename_list):
        command = 'rm %s' % ' '.join(filename for filename in filename_list)
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def clear_actual():
        logging.info("clearing actual...")
        command = 'kdg actual'
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    @staticmethod
    def clear_mos():
        logging.info("clearing mos...")
        command = 'kdg scfmo'
        subprocess.call(command, shell=True)
        command = 'adg scfmo none file=mos'
        subprocess.call(command, shell=True)

    @staticmethod
    def increase_scfiterlimit():
        logging.info("increasing scfiterlimit...")
        command = 'kdg scfiterlimit'
        subprocess.call(command, shell=True)
        command = 'adg scfiterlimit 3000'
        subprocess.call(command, shell=True)

    def read_result(self, energy_type=None, orbital_to_read=None, supplied_folder_path=None):
        """Reads either the resulting total or orbital energy from dscf/ridft output."""

        if supplied_folder_path:
            calc_folder_path = supplied_folder_path
        else:
            calc_folder_path = self.calc_folder_path
        out_file_path = os.path.join(calc_folder_path, '%s.log' % self.calculation_type)
        out_file = open(out_file_path, 'r')
        output_energies = []
        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if energy_type == 'orbital':
                    if split_line[0] == "irrep" and orbital_to_read['irrep'] in split_line:
                        energy_line = ' '.join(linecache.getline(out_file_path, i + 3).split()).split()
                        orbital_e = energy_line[split_line.index(orbital_to_read['irrep'])]
                        logging.info("    new orbital energy: %s" % orbital_e)
                        output_energies.append(orbital_e)
                        linecache.clearcache()
                elif energy_type == 'total':
                    if split_line[0] == "|" and "total" in split_line and "energy" in split_line:
                        energy_line = ' '.join(linecache.getline(out_file_path, i + 1).split()).split()
                        # Correct for the stupid output format that sometimes doesn't leave the spaces around the energy
                        if energy_line[4] == '|':
                            total_e = energy_line[3][1:]
                        else:
                            total_e = energy_line[4]
                        logging.info("    new total energy: %s" % total_e)
                        output_energies.append(total_e)
                        linecache.clearcache()
                else:
                    logging.info('Please specify a valid energy type (orbital or total) for me to optimise.')
        logging.info('   Found energy for type: %s, irrep: %s: %s' % (
            energy_type,
            orbital_to_read['irrep'] if orbital_to_read else 'total',
            float(output_energies[self.spin_indices[orbital_to_read['spin']]]) if orbital_to_read
            else output_energies[0]))
        return float(output_energies[self.spin_indices[orbital_to_read['spin']]]) if orbital_to_read \
            else float(output_energies[0])

    def read_eiger(self, value_to_read='homolumo'):
        self.basis.run_eiger(add_to_log=True, file_path=self.calc_folder_path, shell=True)
        out_file_path = os.path.join(self.calc_folder_path, 'eiger.log')
        out_file = open(out_file_path, 'r')
        gap_energy = None

        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if value_to_read == 'homolumo':
                    if split_line[0] == 'Gap':
                        gap_line = ' '.join(linecache.getline(out_file_path, i + 1).split()).split()
                        gap_energy = gap_line[5]
                        logging.info("    new HOMO-LUMO gap: %s eV" % gap_energy)
                        linecache.clearcache()

        return float(gap_energy)

    def read_escf(self, excitation_to_read=None):
        logging.info('reading excitation: %s...' % excitation_to_read['repr'])
        out_file_path = os.path.join(self.calc_folder_path, '%s.log' % 'escf')
        out_file = open(out_file_path, 'r')
        excitation_energy = 100.0
        oscillator_strength = 0.0
        rotatory_strength = 0.0
        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if split_line[-1] == "excitation" and split_line[0]+split_line[2] == excitation_to_read['repr']:
                    if not 'WARNING!' in linecache.getline(out_file_path, i + 7):
                        logging.info(' ... found excitation')
                        energy_ev = ' '.join(linecache.getline(out_file_path, i + 8).split()).split()[4]
                        oscillator_str = ' '.join(linecache.getline(out_file_path, i + 17).split()).split()[2]
                        rotatory_str = ' '.join(linecache.getline(out_file_path, i + 26).split()).split()[2]

                        logging.info("    new excitation energy: %s eV\n"
                                     "      oscillator strength: %s\n"
                                     "        rotatory strength: %s" % (energy_ev, oscillator_str, rotatory_str))

                        excitation_energy = float(energy_ev)
                        oscillator_strength = float(oscillator_str)
                        rotatory_strength = float(rotatory_str)
                        linecache.clearcache()

        return {'excitation': excitation_energy,
                'oscillation': oscillator_strength,
                'rotation': rotatory_strength
                }

    def run_multivariate_orbital_optimisation(self, array_of_potentials):
        """Tries to optimise for tracked orbital energies using specified potentials."""

        def run_multivariate_calc(x0):
            """Takes array of coeffs and exps, maps them to and updates their ecps, runs the calculation and
               reads the results from calc.log. It then normalises and returns the results as a list."""
            self.iterations += 1
            logging.info('\n============ Optimisation iteration: %s (Seed %s) ============' % (self.iterations, self.current_seed))

            potential_variable_list = []
            if self.optimise_pseudo_geometry:
                for pot in range(0, len(x0)-2, 2):
                    potential_variable_list.append([x0[pot], x0[pot+1]])
            else:
                for pot in range(0, len(x0), 2):
                    potential_variable_list.append([x0[pot], x0[pot+1]])

            # Use list of ecp locators to insert all variables
            if len(potential_variable_list) == len(self.ecp_locators):
                for i, ecp_locator in enumerate(self.ecp_locators):
                    #ecp_locator['functions_list'][0]['coefficient'] = potential_variable_list[i][0]
                    #ecp_locator['functions_list'][0]['exponent'] = potential_variable_list[i][1]
                    ecp_locator['functions_list'] = [{'coefficient': potential_variable_list[i][0],
                                                      'r^n': 1,
                                                      'exponent': potential_variable_list[i][1]}]

                    if self.optimise_with_total_gaps is True:
                        for tracked_total_gap in self.tracked_total_gaps:
                            self.basis.variable_file_path = os.path.join(tracked_total_gap['gap_folder_path'], 'basis')
                            self.basis.update_variable(self.ecp_locators[i])
                            logging.info('Updated ECP %s: coeff %s exp %s' % (ecp_locator['basis_ecp_name'],
                                                                          potential_variable_list[i][0],
                                                                          potential_variable_list[i][1]))
                    if self.calc_folder_path:
                        self.basis.variable_file_path = os.path.join(self.calc_folder_path, 'basis')
                    self.basis.update_variable(self.ecp_locators[i])
                    logging.info('Updated ECP %s: coeff %s exp %s' % (ecp_locator['basis_ecp_name'],
                                                                      potential_variable_list[i][0],
                                                                      potential_variable_list[i][1]))
            else:
                logging.info('Error mapping variables to ecps. No of ecps: %s. No of variable pairs: %s' %
                      (len(self.ecp_locators), len(potential_variable_list)))

            # Update potential geometry
            if self.optimise_pseudo_geometry:
                if self.optimise_with_total_gaps:
                    for tracked_total_gap in self.tracked_total_gaps:
                        self.coord.coord_file_path = os.path.join(tracked_total_gap['gap_folder_path'], 'coord')
                        self.coord.read_coords()
                        for pseudo_carbon in self.pseudo_geometry['indices_of_pseudo_carbons']:
                            if self.pseudo_geometry_type == 'sp2':
                                self.coord.repseudopotentialise_sp2_atom(pseudo_carbon, x0[-2], x0[-1])
                            elif self.pseudo_geometry_type == 'sp3':
                                self.coord.set_potential_distance_to(pseudo_carbon, x0[-2])
                            logging.info('Re-potentialised atom %s, with set distance %s, split %s' % (
                                self.pseudo_geometry['indices_of_pseudo_carbons'],
                                x0[-2],
                                x0[-1]))

                # self.coord.coord_list = []
                self.coord.coord_file_path = os.path.join(self.calc_folder_path, 'coord')
                self.coord.read_coords()
                for pseudo_carbon in self.pseudo_geometry['indices_of_pseudo_carbons']:
                    if self.pseudo_geometry_type == 'sp2':
                        self.coord.repseudopotentialise_sp2_atom(pseudo_carbon, x0[-2], x0[-1])
                    elif self.pseudo_geometry_type == 'sp3':
                        self.coord.set_potential_distance_to(pseudo_carbon, x0[-2])
                    logging.info('Re-potentialised atom %s, with set distance %s, split %s' % (
                        self.pseudo_geometry['indices_of_pseudo_carbons'],
                        x0[-2],
                        x0[-1]))

            self.clear_actual()
            self.clear_mos()
            self.increase_scfiterlimit()
            self.remove_files(['moments'])

            # Run the calculation
            if self.calculation_type == 'dscf':
                self.basis.run_dscf(add_to_log=True, file_path=self.calc_folder_path, shell=True)
            elif self.calculation_type == 'ridft':
                self.basis.run_ridft(add_to_log=True, file_path=self.calc_folder_path, shell=True)

            # Read all results

            # Read orbital energies
            read_orbital_energies = []
            if self.optimise_with_orbitals is True:
                for tracked_orbital in self.tracked_orbitals:
                    read_orbital_energies.append(self.read_result(energy_type='orbital',
                                                                  orbital_to_read=tracked_orbital))
            # Normalise orbitals
            normalised_orbital_energies = []
            if self.optimise_with_orbitals is True:
                for i, read_orbital_energy in enumerate(read_orbital_energies):
                    normalised_energy = numpy.abs(read_orbital_energy - self.tracked_orbitals[i]['reference_energy'])\
                                                     * self.tracked_orbitals[i]['orbital_error_multiplier']
                    logging.info('             irrep %s: %s eV (error %s, multiplier %s)' % (
                        self.tracked_orbitals[i]['irrep'],
                        read_orbital_energy,
                        normalised_energy,
                        self.tracked_orbitals[i]['orbital_error_multiplier']
                        )
                    )
                    normalised_orbital_energies.append(normalised_energy)

            # Normalise total
            normalised_total = 0.0
            if self.optimise_with_total_energy is True:
                read_total_energy = self.read_result(energy_type='total')
                normalised_total = numpy.abs(read_total_energy - self.tracked_total_energy)

            normalised_total_gap_error = 0.0
            if self.optimise_with_total_gaps is True:
                cumulative_total_gap_error = 0.0
                read_this_energy = self.read_result(energy_type='total')
                for tracked_total_gap in self.tracked_total_gaps:
                    # Run the other calc
                    if self.calculation_type == 'dscf':
                        self.basis.run_dscf(add_to_log=True, file_path=tracked_total_gap['gap_folder_path'], shell=True)
                    elif self.calculation_type == 'ridft':
                        self.basis.run_ridft(add_to_log=True, file_path=tracked_total_gap['gap_folder_path'], shell=True)
                    # read the other energy
                    read_that_energy = self.read_result(energy_type='total', supplied_folder_path=tracked_total_gap['gap_folder_path'])
                    # find difference
                    total_gap_energy_eh = read_that_energy - read_this_energy
                    total_gap_energy_ev = total_gap_energy_eh * 27.2113845  # Convert to eV as everything else is in eV
                    total_gap_error = numpy.abs(total_gap_energy_ev - tracked_total_gap['reference_gap_energy'])
                    cumulative_total_gap_error += total_gap_error
                    logging.info('  Total gap between %s and %s: %s eV \n'
                                 '                              (error %s eV)' %
                                 (self.calc_folder_path,
                                  tracked_total_gap['gap_folder_path'],
                                  total_gap_energy_ev,
                                  total_gap_error))
                normalised_total_gap_error = cumulative_total_gap_error * self.total_gaps_error_multiplier
                logging.info('Total total gap error %s eV, multiplier %s' %
                             (normalised_total_gap_error,
                              self.total_gaps_error_multiplier))


            # Collect ESCF results
            if self.optimise_with_fixed_excitations is True \
                    or self.optimise_with_flexible_excitations is True \
                    or self.optimise_with_excitation_spectra is True:
                self.remove_files(['sing*', 'vfile*', 'wfile*'])
                self.basis.run_escf(add_to_log=True, file_path=self.calc_folder_path)

            normalised_excitation_energy = 0.0
            normalised_oscillator_strength = 0.0
            normalised_rotatory_strength = 0.0
            read_escf_results = []
            if self.optimise_with_fixed_excitations is True:
                for tracked_excitation in self.tracked_excitations:
                    read_escf_results.append(self.read_escf(excitation_to_read=tracked_excitation))

                for i, read_excitation in enumerate(read_escf_results):
                    normalised_excitation = numpy.abs(read_excitation['excitation'] - self.tracked_excitations[i]['reference_energy'])
                    normalised_oscillation = numpy.abs(read_excitation['oscillation'] - self.tracked_excitations[i]['reference_oscillation'])
                    normalised_rotation = numpy.abs(read_excitation['rotation'] - self.tracked_excitations[i]['reference_rotation'])
                    logging.info("     new excitation error: %s eV" % normalised_excitation)
                    logging.info("    new oscillation error: %s eV" % normalised_oscillation)
                    logging.info("       new rotation error: %s eV" % normalised_rotation)
                    normalised_excitation_energy += normalised_excitation
                    normalised_oscillator_strength += normalised_oscillation
                    normalised_rotatory_strength += normalised_rotation

            normalised_flexible_excitation_error = 0.0
            if self.optimise_with_flexible_excitations is True:
                escf_control = ESCFControl()
                escf_control.read_escf()
                for tracked_excitation in self.tracked_excitations:
                    normalised_contribution_error, normalised_oscillator_error, calced_excitation = escf_control.evaluate_peak_error(
                        tracked_excitation['reference_energy'],
                        tracked_excitation['reference_oscillation'],
                        tracked_excitation['electric_dipole_norm'],
                        tracked_excitation['repr'])
                    if calced_excitation:
                        logging.info("identified flexible excitation: %s with dipole norm %s" %
                                     (calced_excitation['id'],
                                      calced_excitation['electric_dipole_norm']))
                    else:
                        logging.info("        no matching excitation found.")
                    logging.info("flexible excitation wavelength error: %s eV" % normalised_contribution_error)
                    #logging.info("flexible excitation oscillator error: %s" % normalised_oscillator_error)
                    normalised_flexible_excitation_error += normalised_contribution_error
                    normalised_flexible_excitation_error += normalised_oscillator_error * self.flexible_excitation_osc_error_multiplier

            normalised_excitation_spectra_error = 0.0
            # normalised_excitation_spectral_hash_error = 0.0
            if self.optimise_with_excitation_spectra is True:
                escf_control = ESCFControl()
                spectral_error = escf_control.evaluate_spectral_error()
                normalised_excitation_spectra_error += spectral_error * self.excitation_spectra_error_multiplier
                # normalised_excitation_spectral_hash_error += spectral_hash_error * self.excitation_spectra_error_multiplier
                logging.info("           new spectral least squares error: %s (multiplier %s)" % (spectral_error, self.excitation_spectra_error_multiplier))
                # logging.info("      new spectral hash error: %s eV" % spectral_hash_error)


            # Read HOMO-LUMO gap
            normalised_homo_lumo_gap = 0.0
            if self.optimise_with_homo_lumo_gap is True:
                read_gap_energy = self.read_eiger(value_to_read='homolumo')
                normalised_homo_lumo_gap = numpy.abs(read_gap_energy - self.tracked_homo_lumo_gap)
                logging.info("      new HOMO-LUMO error: %s eV" % normalised_homo_lumo_gap)

            # Read energy gradients
            gradient_error = 0.0
            if len(self.tracked_atom_gradients) > 0:
                self.gradient.run_grad(file_path=self.calc_folder_path)
                self.gradient.variable_file_path = self.calc_folder_path + 'gradient'
                self.gradient.read_variables()
                watched_atom_indices = (index for (index, d) in enumerate(self.gradient.gradient_file[-1]['atoms']) if d['#'] in self.tracked_atom_gradients)
                for atom_index in watched_atom_indices:
                    dx = self.gradient.gradient_file[-1]['atoms'][atom_index]['dx']
                    dy = self.gradient.gradient_file[-1]['atoms'][atom_index]['dy']
                    dz = self.gradient.gradient_file[-1]['atoms'][atom_index]['dz']
                    gradient_error = gradient_error + numpy.abs(dx) + numpy.abs(dy) + numpy.abs(dz)
                gradient_error = gradient_error * self.gradient_error_multiplier
                logging.info('Gradient error: %s (multiplier %s)' % (gradient_error, self.gradient_error_multiplier))

            total_error = sum(normalised_orbital_energies) * self.orbital_error_multiplier \
                + normalised_homo_lumo_gap * self.homo_lumo_gap_error_multiplier \
                + normalised_total_gap_error \
                + gradient_error \
                + normalised_total \
                + normalised_excitation_energy \
                + normalised_flexible_excitation_error * self.flexible_excitation_error_multiplier \
                + normalised_excitation_spectra_error \
                + normalised_oscillator_strength \
                + normalised_rotatory_strength


            logging.info('\nTotal error: %s eV\n' % str(total_error))

            # Plot spectral optimisation
            if self.optimise_with_excitation_spectra:
                escf_control = ESCFControl()
                reference_convolutions = escf_control.read_spyctrum(escf_control.reference_convolution_path)
                pseudo_convolutions = escf_control.read_spyctrum(escf_control.convolution_path)
                x_ref = [float(convolution[0]) for convolution in reference_convolutions]
                x_ps = [float(convolution[0]) for convolution in pseudo_convolutions]
                y_ref = [float(convolution[1]) for convolution in reference_convolutions]
                y_ps = [float(convolution[1]) for convolution in pseudo_convolutions]

                plt.plot(x_ref, y_ref)
                plt.plot(x_ps, y_ps, linestyle="dashed")
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('Intensity of absorption (arb. units)')
                plt.text(900, 2.5, 'Iteration: %s' % "{:05}".format(self.iterations))
                leading_iteration = "{:05}".format(self.iterations)
                plt.savefig('%s/spectracomp_%s.png' % (self.moo_data_directory, leading_iteration))
                plt.gcf().clear()

            return total_error

        logging.info('Optimising potentials with: '
                     '\n  - ecps %s %s %s %s %s %s %s %s %s' % (
                     ', '.join(locator['basis_ecp_name'] for locator in self.ecp_locators),
                     '\n  - over irreps: %s' % ', '.join(tracked['irrep'] for tracked in self.tracked_orbitals) if self.tracked_orbitals else '',
                     '\n  - fixed excitations %s' % ', '.join(tracked['repr'] for tracked in self.tracked_excitations) if self.optimise_with_fixed_excitations else '',
                     '\n  - flexible excitations with dipole norms: %s' % str([tracked['electric_dipole_norm'] for tracked in self.tracked_excitations]) if self.optimise_with_flexible_excitations else '',
                     '\n  - excitation spectra' if self.optimise_with_excitation_spectra else '',
                     '\n  - total energy' if self.optimise_with_total_energy else '',
                     '\n  - HOMO-LUMO gap' if self.optimise_with_homo_lumo_gap else '',
                     '\n  - a total energy gap between this folder and folders %s' % str([tracked['gap_folder_path'] for tracked in self.tracked_total_gaps]) if self.optimise_with_total_gaps else '',
                     '\n  - variable pseudo-potential geometry' if self.optimise_pseudo_geometry else '')
                     )

        subprocess.call('rm -rf %s; mkdir %s' % (self.moo_data_directory, self.moo_data_directory), shell=True)

        #                         coefficient   exponent
        single_potential_bounds = (-50, 50), (0.001, 50)
        all_bounds = []
        for ecp_to_optimise in self.ecp_locators:
            all_bounds.extend(single_potential_bounds)

        # Only SLSQP handles constraints and bounds.

        if self.optimise_pseudo_geometry:
            all_bounds.append(
                (self.pseudo_geometry['potential_set_distance_bounds'][0],
                 self.pseudo_geometry['potential_set_distance_bounds'][1])
            )
            all_bounds.append(
                (self.pseudo_geometry['potential_set_split_distance_bounds'][0],
                 self.pseudo_geometry['potential_set_split_distance_bounds'][1])
            )
            array_of_potentials = numpy.append(array_of_potentials, self.pseudo_geometry['potential_set_distance_guess'])
            array_of_potentials = numpy.append(array_of_potentials, self.pseudo_geometry['potential_set_split_distance_guess'])

        logging.info('Commencing optimisation...')
        try:
            self.coord.read_coords()
            if self.seeded_optimisation:
                logging.info('Running with %s random starting seed(s)... \n'
                             'Mu-weighted variables: %s ' % (self.initial_seed_number, array_of_potentials))
                results_dicts_all_seeds = []

                def create_bounded_weighted_seed(bounds, mu=0):
                    """Most successful potentials are rarely at the bound extremes.
                    Weighting the seeds should cut down on number required."""
                    #mu, sigma = 0, 0.1  # mean and standard deviation
                    #s = np.random.normal(mu, sigma, 1000)
                    while True:
                        seed_guess = numpy.random.normal(mu, 0.5*bounds[1])
                        if bounds[0] < seed_guess < bounds[1]:
                            return seed_guess

                for seed in range(self.initial_seed_number):
                    self.current_seed += 1
                    logging.info('\n============ Seed: %s ============' % self.current_seed)
                    initial_guesses = []
                    for index, bounds_set in enumerate(all_bounds):
                        #variable_guess = numpy.random.uniform(low=bounds_set[0], high=bounds_set[1])
                        variable_guess = create_bounded_weighted_seed(bounds_set, array_of_potentials[index])
                        initial_guesses.append(variable_guess)

                    min_test = scipy.optimize.minimize(run_multivariate_calc,
                                                       initial_guesses,
                                                       bounds=all_bounds,
                                                       # method='Nelder-Mead',
                                                       # method='Powell',
                                                       # method='CG',
                                                       # method='BFGS',
                                                       method='SLSQP',
                                                       # method='COBYLA',
                                                       options={'eps': 0.001,
                                                                'maxiter': 150},
                                                       tol=0.0001
                                                       )

                    min_test['moo_iteration'] = self.iterations
                    min_test['seed_number'] = self.current_seed
                    results_dicts_all_seeds.append(min_test)
                    logging.info(min_test)

                    sorted_results = sorted(results_dicts_all_seeds, key=lambda k: k['fun'])
                    logging.info('Ran optimisation with %s initial seed(s), best result: \n %s' %
                                 (self.current_seed, sorted_results[0]))

                    with open('seed_results.moo', 'w') as seed_result_file:
                        try:
                            seed_result_file.write('%s MOO seed results, smallest to largest total error:\n'
                                                   % self.current_seed)
                            seed_result_file.write('%s \n' % [sorted_result for sorted_result in sorted_results])
                        except Exception as e:
                            seed_result_file.write('%s' % e)
                            seed_result_file.close()

                    leading_iteration = "{:05}".format(self.iterations)
                    command = 'mkdir %s/iteration_%s; cp * %s/iteration_%s' % (self.moo_data_directory,
                                                                               leading_iteration,
                                                                               self.moo_data_directory,
                                                                               leading_iteration)
                    subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            else:
                logging.info('Running with specified starting guesses: ' % array_of_potentials)
                output = scipy.optimize.minimize(run_multivariate_calc, array_of_potentials,
                                                 options={'eps': 0.001,
                                                          'maxiter': 150},
                                                 tol=0.0001,
                                                 bounds=all_bounds,
                                                 # method='Nelder-Mead',
                                                 # method='Powell',
                                                 # method='CG',
                                                 # method='BFGS',
                                                 method='SLSQP',
                                                 # method='COBYLA',
                                                 )
                logging.info(output)

            if self.optimise_with_excitation_spectra:
                import imageio
                filenames = []
                if self.iterations > 10:
                    for i in range(1, self.iterations):
                        filenames.append('%s/spectracomp_' % self.moo_data_directory + "{:05}".format(i) + '.png')
                    images = []
                    for filename in filenames:
                        images.append(imageio.imread(filename))
                    imageio.mimsave('spectraopt.gif', images)

        except Exception as e:
            logging.exception(str(e), exc_info=True)


if __name__ == "__main__":
    optimiser = Optimiser()
    if '-set' in sys.argv:
        optimiser.create_empty_settings_file()
    else:
        optimiser.run()


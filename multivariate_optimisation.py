#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
from chem_lib.gradient import GradientControl
from chem_lib.escf import ESCFControl
import scipy.optimize
import numpy
import json
import subprocess
import sys
import logging

empty_setup_file = {
    'calc_folder_path': '',
    'optimise_with_orbitals': True,
    'optimise_with_total_energy': False,
    'optimise_with_fixed_excitations': False,
    'optimise_with_flexible_excitations': False,
    'optimise_with_homo_lumo_gap': False,
    'ecp_locators': [
        {'line_type': '$ecp',
         'basis_ecp_name': 'c ecp-p',
         'orbital_descriptor': 'p-f',
         'functions_list': [{'coefficient': 0,
                             'r^n': 1,
                             'exponent': 0}]
         },
        {'line_type': '$ecp',
         'basis_ecp_name': 'he ecp-s',
         'orbital_descriptor': 's-f',
         'functions_list': [{'coefficient': 0,
                             'r^n': 1,
                             'exponent': 0}]
         }
    ],
    'tracked_orbitals': [{
        'irrep': '1a',
        'spin': 'mos',
        'reference_energy': 0.0
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
    'initial_guesses': [0.1, 0.1, 0.1, 0.1]
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
        self.basis = BasisControl()
        self.gradient = GradientControl()
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.calculation_type = 'dscf'
        self.calc_folder_path = None
        self.iterations = 0
        # optimisation specs
        self.ecp_locators = []

        self.optimise_with_total_energy = False
        self.tracked_total_energy = 0.0

        self.optimise_with_orbitals = False
        self.tracked_orbitals = []  # locators for orbitals to attempt to optimise
        self.orbital_error_multiplier = 1.0

        self.optimise_with_fixed_excitations = False
        self.optimise_with_flexible_excitations = False
        self.tracked_excitations = []  # locators for tddft excitations to attempt to optimise
        self.flexible_excitation_error_multiplier = 1.0
        self.flexible_excitation_osc_error_multiplier = 1.0

        self.tracked_atom_gradients = []  # hashes of atoms from which to read gradient
        self.gradient_error_multiplier = 10.0  # multiplier gradient minimisation

        self.optimise_with_homo_lumo_gap = False
        self.tracked_homo_lumo_gap = 0.0
        self.homo_lumo_gap_error_multiplier = 1.0

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
                                              //    //  
                                              ^^    ^^  
            #######################################
            ''')
        opt_data_path = os.path.join(os.getcwd(), self.setup_file_name)
        if not os.path.isfile(opt_data_path):
            logging.info('No moo file, running setup...')
            self.setup_menu()
        else:
            logging.info('moo file found.')
        with open(opt_data_path, 'r') as opt_data_file:
            optdata = json.load(opt_data_file)
        opt_data_file.close()
        self.optimise_with_orbitals = optdata['optimise_with_orbitals']
        self.optimise_with_total_energy = optdata['optimise_with_total_energy']
        self.optimise_with_fixed_excitations = optdata['optimise_with_fixed_excitations']
        self.optimise_with_flexible_excitations = optdata['optimise_with_flexible_excitations']
        self.optimise_with_homo_lumo_gap = optdata['optimise_with_homo_lumo_gap']
        self.tracked_homo_lumo_gap = optdata['tracked_homo_lumo_gap']
        self.tracked_orbitals = optdata['tracked_orbitals']
        self.tracked_excitations = optdata['tracked_excitations']
        self.tracked_total_energy = optdata['tracked_total']
        self.calc_folder_path = optdata['calc_folder_path']
        self.ecp_locators = optdata['ecp_locators']
        self.initial_guesses = numpy.array(optdata['initial_guesses'])
        if '-nosetup' in sys.argv:
            optimised_value = self.run_multivariate_orbital_optimisation(numpy.array(self.initial_guesses))
            logging.info(optimised_value)
        else:
            if input('Run now? y/n: ') == 'y':
                optimised_value = self.run_multivariate_orbital_optimisation(numpy.array(self.initial_guesses))
                logging.info(optimised_value)

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
        print("opt.moo settings file created. Use this file to add optimisation"
              "of excitations, HOM-LUMO gap, total energy or gradients.")

    def remove_files(self, filename_list):
        command = 'rm %s' % ' '.join(filename for filename in filename_list)
        subprocess.call(command, shell=True)

    @staticmethod
    def clear_actual():
        logging.info("clearing actual...")
        command = 'kdg actual'
        subprocess.call(command, shell=True)

    @staticmethod
    def clear_mos():
        logging.info("clearing mos...")
        command = 'kdg scfmo'
        subprocess.call(command, shell=True)
        command = 'adg scfmo none file=mos'
        subprocess.call(command, shell=True)

    def read_result(self, energy_type=None, orbital_to_read=None):
        """Reads either the resulting total or orbital energy from dscf/ridft output."""

        out_file_path = os.path.join(self.calc_folder_path, '%s.log' % self.calculation_type)
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
            logging.info('\n===== Optimisation iteration: %s =====' % self.iterations)

            potential_variable_list = []
            for pot in range(0, len(x0), 2):
                potential_variable_list.append([x0[pot], x0[pot+1]])

            # Use list of ecp locators to insert all variables
            if len(potential_variable_list) == len(self.ecp_locators):
                for i, ecp_locator in enumerate(self.ecp_locators):
                    ecp_locator['functions_list'][0]['coefficient'] = potential_variable_list[i][0]
                    ecp_locator['functions_list'][0]['exponent'] = potential_variable_list[i][1]
                    if self.calc_folder_path:
                        self.basis.variable_file_path = os.path.join(self.calc_folder_path, 'basis')
                    self.basis.update_variable(self.ecp_locators[i])
                    logging.info('Updated ECP %s: coeff %s exp %s' % (ecp_locator['basis_ecp_name'],
                                                                      potential_variable_list[i][0],
                                                                      potential_variable_list[i][1]))
            else:
                logging.info('Error mapping variables to ecps. No of ecps: %s. No of variable pairs: %s' %
                      (len(self.ecp_locators), len(potential_variable_list)))

            self.clear_actual()
            self.clear_mos()

            # Run the calculation
            if self.calculation_type == 'dscf':
                self.basis.run_dscf(add_to_log=True, file_path=self.calc_folder_path, shell=True)
            elif self.calculation_type == 'ridft':
                self.basis.run_ridft(add_to_log=True, file_path=self.calc_folder_path, shell=True)

            # Read all results

            # Read orbital energies
            read_orbital_energies = []
            for tracked_orbital in self.tracked_orbitals:
                read_orbital_energies.append(self.read_result(energy_type='orbital',
                                                              orbital_to_read=tracked_orbital))
            # Normalise orbitals
            normalised_orbital_energies = []
            if self.optimise_with_orbitals is True:
                for i, read_orbital_energy in enumerate(read_orbital_energies):
                    normalised_energy = numpy.abs(read_orbital_energy - self.tracked_orbitals[i]['reference_energy'])
                    logging.info('             irrep %s: %s eV (error %s)' % (
                        self.tracked_orbitals[i]['irrep'],
                        read_orbital_energy,
                        normalised_energy
                        )
                    )
                    normalised_orbital_energies.append(normalised_energy)

            # Normalise total
            normalised_total = 0.0
            if self.optimise_with_total_energy is True:
                read_total_energy = self.read_result(energy_type='total')
                normalised_total = numpy.abs(read_total_energy - self.tracked_total_energy)

            # Collect ESCF results
            if self.optimise_with_fixed_excitations is True or self.optimise_with_flexible_excitations is True:
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
                    #normalised_flexible_excitation_error += normalised_oscillator_error * self.flexible_excitation_osc_error_multiplier

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
                + gradient_error \
                + normalised_total \
                + normalised_excitation_energy \
                + normalised_flexible_excitation_error * self.flexible_excitation_error_multiplier \
                + normalised_oscillator_strength \
                + normalised_rotatory_strength


            logging.info('\nTotal error: %s eV\n' % str(total_error))

            command = 'mkdir moo_iterations/iteration_%s; cp * moo_iterations/iteration_%s' % (self.iterations, self.iterations)
            subprocess.call(command, shell=True)

            return total_error

        logging.info('Optimising potentials with: '
                     '\n  - ecps %s over irreps %s %s %s %s %s' % (
                     ', '.join(locator['basis_ecp_name'] for locator in self.ecp_locators),
                     ', '.join(tracked['irrep'] for tracked in self.tracked_orbitals),
                     '\n  - excitations %s' % ', '.join(tracked['repr'] for tracked in self.tracked_excitations) if self.optimise_with_fixed_excitations else '',
                     '\n  - flexible excitations with dipole norms: %s' % str([tracked['electric_dipole_norm'] for tracked in self.tracked_excitations]) if self.optimise_with_flexible_excitations else '',
                     '\n  - total energy' if self.optimise_with_total_energy else '',
                     '\n  - HOMO-LUMO gap' if self.optimise_with_homo_lumo_gap else '')
        )

        subprocess.call('rm -rf moo_iterations; mkdir moo_iterations', shell=True)

        # Bounds:   p coeff -ve,   p exp +ve,   s coeff +ve,   s exp +ve
        pp_bounds = [(-50, 50), (0.001, 50), (-50, 50), (0.001, 50)]
        #pp_bounds = [(-50, 50), (0.001, 50), (-50, 50), (0.001, 50), (-100, 100), (0.001, 50)]
        # Only SLSQP handles constraints and bounds.
        logging.info('Commencing optimsation...')
        try:
            return scipy.optimize.minimize(run_multivariate_calc, array_of_potentials,
                                           options={'eps': 0.001,
                                                    'maxiter': 150},
                                           tol=0.00001,
                                           bounds=pp_bounds,
                                           # method='Nelder-Mead',
                                           # method='Powell',
                                           # method='CG',
                                           # method='BFGS',
                                           method='SLSQP',
                                           # method='COBYLA',
                                           )
        except Exception as e:
            logging.exception(str(e), exc_info=True)


if __name__ == "__main__":
    optimiser = Optimiser()
    optimiser.run()


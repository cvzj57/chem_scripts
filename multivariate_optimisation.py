#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
from chem_lib.gradient import GradientControl
import scipy.optimize
import numpy
import json

empty_setup_file = {
    'calc_folder_path': '',
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
    'tracked_orbitals': [
        {'irrep': '1a'},
        {'spin': 'mos'},
        {'reference_energy': 0.0},
    ],
    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)
    'initial_guesses': [0.1, 0.1, 0.1, 0.1]
}


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energy."""
    def __init__(self):
        self.reference_E = None
        self.setup_file_name = 'opt.moo'
        self.basis = BasisControl()
        self.gradient = GradientControl()
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.calculation_type = 'dscf'
        self.calc_folder_path = None
        self.tracked_orbitals = []  # locators for orbitals to attempt to optimise
        self.tracked_atom_gradients = []  # hashes of atoms from which to read gradient
        self.ecp_locators = []  # initial guesses for optimiser
        self.gradient_error_multiplier = 10.0  # multiplier gradient minimisation
        self.initial_guesses = None
        with open(os.path.dirname(os.path.realpath(__file__))+'/chem_lib/optimiser_orbital_library.json', 'r') as json_file:
            self.orbital_library = json.load(json_file)

    def run(self):
        print(
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
            print('No moo file, running setup...')
            self.setup_menu()
        else:
            print('moo file found.')
        with open(opt_data_path, 'r') as opt_data_file:
            optdata = json.load(opt_data_file)
        opt_data_file.close()
        self.tracked_orbitals = optdata['tracked_orbitals']
        self.calc_folder_path = optdata['calc_folder_path']
        self.ecp_locators = optdata['ecp_locators']
        self.initial_guesses = numpy.array(optdata['initial_guesses'])
        if input('Run now? y/n: ') == 'y':
            optimised_value = self.run_multivariate_orbital_optimisation(numpy.array(self.initial_guesses))
            print(optimised_value)

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
            setup_file['ecp_locators'].append(
                {'line_type': '$ecp',
                 'basis_ecp_name': 'c ecp-p',
                 'orbital_descriptor': 's-f',
                 'functions_list': [{'coefficient': 0,
                                     'r^n': 1,
                                     'exponent': 0}]
                 })

        with open(os.path.join(setup_file['calc_folder_path'], self.setup_file_name), 'w') as outfile:
            json.dump(setup_file, outfile, sort_keys=True, indent=2)
        outfile.close()
        print("Whaddaya know, you've got potential.")

    def read_result(self, energy_type=None, orbital_to_read=None):
        """Updates basis file with new coeff/exp value, runs the calculation and reads either the
        resulting total or orbital energy."""

        out_file_path = os.path.join(self.calc_folder_path, 'calc.log')
        out_file = open(out_file_path, 'r')
        output_energies = []
        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if energy_type == 'orbital':
                    if split_line[0] == "irrep" and orbital_to_read['irrep'] in split_line:
                        energy_line = ' '.join(linecache.getline(out_file_path, i + 3).split()).split()
                        orbital_e = energy_line[split_line.index(orbital_to_read['irrep'])]
                        print("new orbital energy %s" % orbital_e)
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
                        print("new total energy %s" % total_e)
                        output_energies.append(total_e)
                        linecache.clearcache()
                else:
                    print('Please specify a valid energy type (orbital or total) for me to optimise.')
        print('Found energy for type: %s, irrep: %s: %s' % (
                       energy_type,
                       orbital_to_read['irrep'],
                       float(output_energies[self.spin_indices[orbital_to_read['spin']]])
                   )
              )
        return float(output_energies[self.spin_indices[orbital_to_read['spin']]])

    def run_multivariate_orbital_optimisation(self, array_of_potentials):
        """Tries to optimise for tracked orbital energies using specified potentials."""

        def run_multivariate_calc(x0):
            '''Takes array of coeffs and exps, maps them to and updates their ecps, runs the calculation and
               reads the results from calc.log. It then normalises and returns the results as a list.'''
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
            else:
                print('Error mapping variables to ecps. No of ecps: %s. No of variable pairs: %s' %
                      (len(self.ecp_locators), len(potential_variable_list)))

            # Run the calculation
            if self.calculation_type == 'dscf':
                self.basis.run_dscf(add_to_log=True, file_path=self.calc_folder_path)
            elif self.calculation_type == 'ridft':
                self.basis.run_ridft(add_to_log=True, file_path=self.calc_folder_path)

            # Read orbital energies
            read_orbital_energies = []
            for tracked_orbital in self.tracked_orbitals:
                read_orbital_energies.append(self.read_result(energy_type='orbital',
                                                              orbital_to_read=tracked_orbital))
            print('Read new energies...')

            # Normalise results
            normalised_orbital_energies = []
            for i, read_orbital_energy in enumerate(read_orbital_energies):
                normalised_energy = numpy.abs(read_orbital_energy - self.tracked_orbitals[i]['reference_energy'])
                print('irrep: %s, %s eV (error %s)' % (
                    self.tracked_orbitals[i]['irrep'],
                    read_orbital_energy,
                    normalised_energy
                    )
                )
                normalised_orbital_energies.append(normalised_energy)

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
                print('Gradient error: %s (multiplier %s)' % (gradient_error, self.gradient_error_multiplier))

            print('Total error: %s eV' % str(sum(normalised_orbital_energies) + gradient_error))
            return sum(normalised_orbital_energies) + gradient_error

        print('Optimising ecps %s over irreps %s' % (
            [locator['basis_ecp_name'] for locator in self.ecp_locators], [tracked['irrep'] for tracked in self.tracked_orbitals])
        )

        # Bounds:   p coeff -ve,   p exp +ve,   s coeff +ve,   s exp +ve
        pp_bounds = [(-50, 50), (0.001, 50), (-50, 50), (0.001, 50)]
        pp_bounds = [(-50, 50), (0.001, 50), (-50, 50), (0.001, 50), (-100, 100), (0.001, 50)]
        # Only SLSQP handles constraints and bounds.
        return scipy.optimize.minimize(run_multivariate_calc, array_of_potentials,
                                       options={'eps': 0.001,
                                                'maxiter': 1000},
                                       tol=0.0000000001,
                                       bounds=pp_bounds,
                                       # method='Nelder-Mead',
                                       # method='Powell',
                                       # method='CG',
                                       # method='BFGS',
                                       method='SLSQP',
                                       # method='COBYLA',
                                       )


if __name__ == "__main__":
    optimiser = Optimiser()
    optimiser.run()


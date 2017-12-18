#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
from chem_lib.gradient import GradientControl
import scipy.optimize
import numpy
import json


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energy."""
    def __init__(self):
        self.reference_E = None
        self.basis = BasisControl()
        self.gradient = GradientControl()
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.calculation_type = 'dscf'
        self.calc_folder_path = 'pseudo_optimisation'
        self.tracked_orbitals = []  # locators for orbitals to attempt to optimise
        self.tracked_atom_gradients = []  # hashes of atoms from which to read gradient
        self.ecp_locators = []
        self.gradient_error_multiplier = 10.0  # multiplier gradient minimisation
        with open('chem_lib/optimiser_orbital_library.json', 'r') as json_file:
            self.orbital_library = json.load(json_file)

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
        pp_bounds = [(None, 0.0), (0.001, None), (None, None), (0.001, None)]
        # pp_bounds = [(None, 0.0), (0.001, None), (0.0, None), (0.001, None), (0.0, None), (0.001, None)]
        # Only SLSQP handles constraints and bounds.
        return scipy.optimize.minimize(run_multivariate_calc, array_of_potentials,
                                       tol=0.0001,
                                       bounds=pp_bounds,
                                       method='Nelder-Mead',
                                       # method='Powell',
                                       # method='CG',
                                       # method='BFGS',
                                       # method='SLSQP',
                                       )


def optimise_energy():

    standard_ecp_locators = [
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
    ]

    methane_3_optimiser = Optimiser()
    ### Methane ecp locators
    methane_3_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    methane_3_optimiser.tracked_orbitals = [
        {'irrep': '1a',
         'spin': 'mos',
         'reference_energy': -14.892},
    ]

    methane_2_optimiser = Optimiser()
    ### Methane ecp locators
    methane_2_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    methane_2_optimiser.tracked_orbitals = [
        {'irrep': '1a1',
         'spin': 'mos',
         'reference_energy': -14.892},
        {'irrep': '1b2',
         'spin': 'mos',
         'reference_energy': -14.892},
    ]

    methane_2_optimiser.tracked_atom_gradients = [1]

    methane_1_optimiser = Optimiser()
    ### Methane ecp locators
    methane_1_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    methane_1_optimiser.tracked_orbitals = [
        {'irrep': '1a',
         'spin': 'mos',
         'reference_energy': -14.892},
        {'irrep': '2a',
         'spin': 'mos',
         'reference_energy': -14.892},
        {'irrep': '3a',
         'spin': 'mos',
         'reference_energy': -14.892},
    ]

    methane_1_lower_optimiser = Optimiser()
    ### Methane ecp locators
    methane_1_lower_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    methane_1_lower_optimiser.tracked_orbitals = [
        {'irrep': '1a',
         'spin': 'mos',
         'reference_energy': -25.350},
        {'irrep': '2a',
         'spin': 'mos',
         'reference_energy': -14.892},
        {'irrep': '3a',
         'spin': 'mos',
         'reference_energy': -14.892},
    ]

    ethene_optimiser = Optimiser()
    ### Ethene ecp locators
    ethene_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    ethene_optimiser.tracked_orbitals = [
        {'irrep': '2a"',
         'spin': 'alpha',
         'reference_energy': -6.630},
        {'irrep': '1a"',
         'spin': 'alpha',
         'reference_energy': -14.512},
    ]

    ethane_optimiser = Optimiser()
    ### Ethane ecp locators
    ethane_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    ethane_optimiser.tracked_orbitals = [
        {'irrep': '5a',
         'spin': 'mos',
         'reference_energy': -14.028},
        {'irrep': '6a',
         'spin': 'mos',
         'reference_energy': -13.328},
        {'irrep': '7a',
         'spin': 'mos',
         'reference_energy': -13.328},
    ]

    ethane_optimiser.tracked_atom_gradients = [1, 2]

    propane_optimiser = Optimiser()
    ### Propane_c ecp locators
    propane_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    propane_optimiser.tracked_orbitals = [
        # {'irrep': '4b2',
        #  'spin': 'mos',
        #  'reference_energy': -9.524},
        # {'irrep': '1a2',
        #  'spin': 'mos',
        #  'reference_energy': -10.609},
        {'irrep': '4a1',
         'spin': 'mos',
         'reference_energy': -9.426},
        {'irrep': '1b1',
         'spin': 'mos',
         'reference_energy': -9.106},
    ]

    half_ethene_optimiser = Optimiser()
    ### half-Ethene ecp locators
    half_ethene_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    half_ethene_optimiser.tracked_orbitals = [
        # {'irrep': '2b2',
        #  'spin': 'mos',
        #  'reference_energy': 7.0},
        {'irrep': '1b2',
         'spin': 'mos',
         'reference_energy': -10.363},
        {'irrep': '1b1',
         'spin': 'mos',
         'reference_energy': -13.771},
        # {'irrep': '3a1',
        #  'spin': 'mos',
        #  'reference_energy': -16.126},
    ]

    #half_ethene_optimiser.tracked_atom_gradients = [1, 2]

    eclipsed_ethane_optimiser = Optimiser()
    ### Eclipsed Ethane ecp locators
    eclipsed_ethane_optimiser.ecp_locators = [
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
         },
        {'line_type': '$basis',
         'basis_ecp_name': 'c def2-svpr',
         'orbital_descriptor': '1s',
         'functions_list': [{'coefficient': 0,
                             'r^n': 1,
                             'exponent': 0}]
         }
    ]

    ## Orbitals to optimise
    eclipsed_ethane_optimiser.tracked_orbitals = [
        {'irrep': '3a1',
         'spin': 'mos',
         'reference_energy': -14.008},
        {'irrep': '1e',
         'spin': 'mos',
         'reference_energy': -13.284},
        {'irrep': '1e',
         'spin': 'mos',
         'reference_energy': -13.284},
    ]

    propane_c_optimiser = Optimiser()
    ### Propane_c ecp locators
    propane_c_optimiser.ecp_locators = standard_ecp_locators

    ## Orbitals to optimise
    propane_c_optimiser.tracked_orbitals = [
        {'irrep': '1b2',
         'spin': 'mos',
         'reference_energy': -14.892},
        {'irrep': '1a1',
         'spin': 'mos',
         'reference_energy': -14.892},
        # {'irrep': '1 b2',
        #  'spin': 'mos',
        #  'reference_energy': -277.846},
    ]


    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)

    # Paola's ethane ecp
    # initial_guesses = numpy.array([-2.8056200103, 0.6244618784275526, 0.5, 1.5])

    # My ethane ecp
    # initial_guesses = numpy.array([-2.8056200103, 0.6, 3.0, 5.5])
    #initial_guesses = numpy.array([-2.090443639, 0.3974474906, 35.58826369, 15.19887795, 0.4, 1.0])
    # methane 2 pot initial
    # initial_guesses = numpy.array([-4.8484, 0.7062, -0.317, 0.3294])
    # initial_guesses = numpy.array([-3.5, 0.8, -0.3, 0.35])
    # initial_guesses = numpy.array([-5.66, 0.88, -0.19, 0.25])
    #initial_guesses = numpy.array([-4.42, 1.86, 0.57, 0.62])

    # Methane 1 pot
    # initial_guesses = numpy.array([-5.25, 2.04, 4.17, 0.43])
    # Methane 1 pot lower
    #initial_guesses = numpy.array([-5.19, 1.81, 0.27, 0.75])

    # basis set 4 (ethene)
    # initial_guesses = numpy.array([-3.9096200103, 0.6244618784, 1.5, 0.5])
    # half pseudo ethene
    initial_guesses = numpy.array([-6.45, 0.90, 1.14, 1.63])
    # initial opt
    # initial_guesses = numpy.array([-8.26, 1.68, -0.95, 0.48])
    # post bopt sopt
    #initial_guesses = numpy.array([-7.2, 1.41, -0.3, 2.1])

    # propane guesses
    # initial_guesses = numpy.array([10.27, 0.29, 100.38, 1.3])

    # base guess
    # initial_guesses = numpy.array([0.1, 0.1, 0.1, 0.1])

    # 4 orbital base guess
    # initial_guesses = numpy.array([-0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    # optimised_value = ethene_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = ethane_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    #optimised_value = eclipsed_ethane_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = propane_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = propane_c_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    optimised_value = half_ethene_optimiser.run_multivariate_orbital_optimisation(initial_guesses)

    # optimised_value = methane_1_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = methane_1_lower_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = methane_2_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = methane_3_optimiser.run_multivariate_orbital_optimisation(initial_guesses)

    print(optimised_value)


optimise_energy()



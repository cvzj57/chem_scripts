#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
import scipy.optimize
import numpy


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energy."""
    def __init__(self):
        self.reference_E = None
        self.basis = BasisControl()
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.calculation_type = 'dscf'
        self.calc_folder_path = 'essai_3'
        self.tracked_orbitals = []
        self.ecp_locators = []

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

            # Read results
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

            # self.normalised_energies
            print('Total error: %s eV' % sum(normalised_orbital_energies))
            return sum(normalised_orbital_energies)

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
    methane_optimiser = Optimiser()
    ### Ethene ecp locators
    methane_optimiser.ecp_locators = [
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
        # {'line_type': '$ecp',
        #  'basis_ecp_name': 'he ecp2-s',
        #  'orbital_descriptor': 's-f',
        #  'functions_list': [{'coefficient': 0,
        #                      'r^n': 1,
        #                      'exponent': 0}]
        #  }
    ]

    ## Orbitals to optimise
    methane_optimiser.tracked_orbitals = [
        {'irrep': '4a',
         'spin': 'mos',
         'reference_energy': -10.923},
        {'irrep': '3a',
         'spin': 'mos',
         'reference_energy': -10.923},
        {'irrep': '2a',
         'spin': 'mos',
         'reference_energy': -19.105},
        # {'irrep': '1a',
        #  'spin': 'mos',
        #  'reference_energy': -277.846},
    ]

    ethene_optimiser = Optimiser()
    ### Ethene ecp locators
    ethene_optimiser.ecp_locators = [
        {'line_type': '$ecp',
         'basis_ecp_name': 'c ecp-p',
         'orbital_descriptor': 'p-f',
         'functions_list': [{'coefficient': 0,
                             'r^n': 1,
                             'exponent': 0}]
         },
        {'line_type': '$ecp',
         'basis_ecp_name': 'h ecp-s',
         'orbital_descriptor': 's-f',
         'functions_list': [{'coefficient': 0,
                             'r^n': 1,
                             'exponent': 0}]
         }
    ]

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
    ethane_optimiser.ecp_locators = [
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

    ## Orbitals to optimise
    ethane_optimiser.tracked_orbitals = [
        {'irrep': '3a',
         'spin': 'mos',
         'reference_energy': -14.028},
        {'irrep': '4a',
         'spin': 'mos',
         'reference_energy': -13.328},
        # {'irrep': '5a',
        #  'spin': 'mos',
        #  'reference_energy': -13.328},
    ]

    propane_c_optimiser = Optimiser()
    ### Propane_c ecp locators
    propane_c_optimiser.ecp_locators = [
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

    ## Orbitals to optimise
    propane_c_optimiser.tracked_orbitals = [
        {'irrep': '1b1',
         'spin': 'mos',
         'reference_energy': -10.923},
        {'irrep': '1a1',
         'spin': 'mos',
         'reference_energy': -19.105},
        # {'irrep': '1 b2',
        #  'spin': 'mos',
        #  'reference_energy': -277.846},
    ]

    # Array of initial guesses (MUST be same order as ecp_locators e.g. p_coeff, p_exp, s_coeff, s_exp)

    # Paola's ethane ecp
    # initial_guesses = numpy.array([-2.8056200103, 0.6244618784275526, 0.5, 1.5])

    # My ethane ecp
    initial_guesses = numpy.array([
        -2.8056200103, 0.6,
        3.0, 5.5
    ])

    # basis set 4 (ethene)
    # initial_guesses = numpy.array([-3.9096200103, 0.6244618784, 1.5, 0.5])

    # base guess
    # initial_guesses = numpy.array([0.1, 0.1, 0.1, 0.1])

    # 4 orbital base guess
    # initial_guesses = numpy.array([-0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    # optimised_value = ethene_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    optimised_value = ethane_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = methane_optimiser.run_multivariate_orbital_optimisation(initial_guesses)
    # optimised_value = propane_c_optimiser.run_multivariate_orbital_optimisation(initial_guesses)

    print(optimised_value)

optimise_energy()



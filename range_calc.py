#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
from scipy.optimize import brentq as brent_opt
from scipy.optimize import brute as brute
# from scipy.optimize import minimize
import scipy
import numpy


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energy."""
    def __init__(self):
        self.reference_E = None
        self.basis = BasisControl()
        self.basis.variable_file_path = 'basis'
        self.trial_range = (-500.0, 1000.0)
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.optimisation_parameter = 'coefficient'
        self.calculation_type = 'dscf'

    def find_energy_for_value(self, trial_value, orbital_to_optimise, energy_type='orbital', cmd_path=''):
        orbital_to_optimise['functions_list'][0][self.optimisation_parameter] = str(trial_value)
        if cmd_path:
            self.basis.variable_file_path = os.path.join(cmd_path, 'basis')
        self.basis.update_variable(orbital_to_optimise)
        if self.calculation_type == 'dscf':
            self.basis.run_dscf(add_to_log=True, file_path=cmd_path)
        elif self.calculation_type == 'ridft':
            self.basis.run_ridft(add_to_log=True, file_path=cmd_path)

        out_file_path = os.path.join(cmd_path, 'calc.log')
        out_file = open(out_file_path, 'r')
        output_energies = []
        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if energy_type == 'orbital':
                    if split_line[0] == "irrep" and orbital_to_optimise['orbital_irrep'] in split_line:
                        energy_line = ' '.join(linecache.getline(out_file_path, i + 3).split()).split()
                        orbital_e = energy_line[split_line.index(orbital_to_optimise['orbital_irrep'])]
                        print("new energy %s eV" % orbital_e)
                        output_energies.append(orbital_e)
                        linecache.clearcache()
                elif energy_type == 'total':
                    if split_line[0] == "|" and "total" in split_line and "energy" in split_line:
                        energy_line = ' '.join(linecache.getline(out_file_path, i + 1).split()).split()
                        total_e = energy_line[4]
                        print("new energy %s eV" % total_e)
                        output_energies.append(total_e)
                        linecache.clearcache()
        print(output_energies)
        return float(output_energies[self.spin_indices[orbital_to_optimise['spin']]])

    def run_minimise_total_e_difference(self,
                                        optimisation_parameter,
                                        orbital_to_optimise,
                                        difference_orbital_to_optimise=None,
                                        reference_energy=None):
        def find_e_difference(trial_value):
            if difference_orbital_to_optimise:
                first_energy = self.find_energy_for_value(trial_value,
                                                          orbital_to_optimise,
                                                          energy_type='total',
                                                          cmd_path='singlet')
                second_energy = self.find_energy_for_value(trial_value,
                                                           difference_orbital_to_optimise,
                                                           energy_type='total',
                                                           cmd_path='triplet')
                difference = first_energy - second_energy
                zeroed = difference - self.reference_E
                print('%s - %s - %s = %s' % (first_energy, second_energy, self.reference_E, zeroed))
                return first_energy - second_energy - self.reference_E
            else:
                return self.find_energy_for_value(trial_value, orbital_to_optimise) - self.reference_E
        self.optimisation_parameter = optimisation_parameter
        self.reference_E = reference_energy
        print('Optimising %s over %s, f(%s) = %s and f(%s) = %s' % (
            optimisation_parameter,
            self.trial_range,
            self.trial_range[0],
            find_e_difference(self.trial_range[0]),
            self.trial_range[1],
            find_e_difference(self.trial_range[1]))
        )
        return scipy.optimize.minimize_scalar(find_e_difference)
        # return scipy.optimize.minimize(find_e_difference, numpy.array([self.trial_range[0]]))
        # return brent_opt(find_e_difference, self.trial_range[0], self.trial_range[1])
        # return brute(find_e_difference, (self.trial_range[0], self.trial_range[1]))


def optimise_energy():
    optimiser = Optimiser()
    orbital_to_optimise = {'line_type': '$ecp',
                           'basis_ecp_name': 'h ecp-s',
                           'orbital_descriptor': 's-f',
                           'orbital_irrep': '1b1u',
                           'spin': 'mos',
                           'functions_list': [{'coefficient': -7.524, 'r^n': 1, 'exponent': 20.0}]
                           }
    difference_orbital_to_optimise = {
                                     'line_type': '$ecp',
                                     'basis_ecp_name': 'h ecp-s',
                                     'orbital_descriptor': 's-f',
                                     'orbital_irrep': '1b1u',
                                     'spin': 'mos',
                                     'functions_list': [{'coefficient': -7.524, 'r^n': 1, 'exponent': 20.0}]
                                     }
    # optimiser.orbital_to_optimise = {'line_type': '$ecp',
    #                                  'basis_ecp_name': 'c ecp-p',
    #                                  'orbital_descriptor': 'p-f',
    #                                  'orbital_irrep': '1a2\"',
    #                                  'spin': 'alpha',
    #                                  'functions_list': [{'coefficient': 3.326, 'r^n': 1, 'exponent': 0.295}]
    #                                  }
    reference_energy = -15.3626756934  # a HF singlet - triplet difference ethene
    # reference_energy = -10.363  # HF ethene pi
    # reference_energy = -6.629  # DFT ethene pi
    print(orbital_to_optimise)

    optimised_value = optimiser.run_minimise_total_e_difference('coefficient',
                                                                orbital_to_optimise,
                                                                difference_orbital_to_optimise=difference_orbital_to_optimise,
                                                                reference_energy=reference_energy)
    print('Optimisation converged to energy %s eV, with %s %s' % (reference_energy,
                                                                  optimiser.optimisation_parameter,
                                                                  optimised_value.x))

optimise_energy()



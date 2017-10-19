#!//share/programs/PYTHON3/bin/python3.6
import linecache
import os
from chem_lib.optimise import BasisControl
from r import MatrixHandler
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
        self.basis.variable_file_path = ''
        self.trial_range = (0.001, 10.0)
        self.spin_indices = {'alpha': -2, 'beta': -1, 'mos': -1}
        self.optimisation_parameter = 'coefficient'
        self.calculation_type = 'dscf'
        self.calc_folder_path = 'try_3_opt'
        self.difference_calc_folder_path = 'triplet'

    def find_energy_for_value(self, trial_value, orbital_to_optimise,
                              energy_type=None, cmd_path=''):
        """Updates basis file with new coeff/exp value, runs the calculation and reads either the
        resulting total or orbital energy."""

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
                       orbital_to_optimise['orbital_irrep'],
                       float(output_energies[self.spin_indices[orbital_to_optimise['spin']]])
                   )
              )
        return float(output_energies[self.spin_indices[orbital_to_optimise['spin']]])

    def run_minimise_total_e_difference(self,
                                        optimisation_parameter,
                                        orbital_to_optimise,
                                        difference_orbital_to_optimise=None,
                                        mos_info=None,
                                        basis_info=None,
                                        reference_energy=None,
                                        use_zeff_calc=False,
                                        energy_type=None):
        """Tries to find root/minimise the energy of a supplied orbital, or the 1st - 2nd
        energy difference between the supplied orbitals. Can be supplied with a trial range and reference energy."""
        def find_e_difference(trial_value, zeff_calc=use_zeff_calc):
            if difference_orbital_to_optimise:

                if zeff_calc:
                    basis_info['file_path'] = os.path.join(self.calc_folder_path, 'basis')
                    mos_info['file_path'] = os.path.join(self.calc_folder_path, 'mos')
                    handler = MatrixHandler()
                    dz_eff = handler.run_handler({'basis_name': 'c def2-SV(P)',
                                                  'orbital_to_extract': 'p'},
                                                 {'orbital': '1 b2',
                                                  'coeff_range': [0, 2]},
                                                 pp_exponent=trial_value,
                                                 basis_path=os.path.join('chem_scripts', 'tm_files', 'basis'),
                                                 mos_path=os.path.join('chem_scripts', 'tm_files', 'alpha'))['dZ_eff']
                    orbital_to_optimise['functions_list'][0]['coefficient'] = -dz_eff
                    difference_orbital_to_optimise['functions_list'][0]['coefficient'] = -dz_eff

                first_energy = self.find_energy_for_value(trial_value,
                                                          orbital_to_optimise,
                                                          energy_type=energy_type,
                                                          cmd_path=self.calc_folder_path)

                second_energy = self.find_energy_for_value(trial_value,
                                                           difference_orbital_to_optimise,
                                                           energy_type=energy_type,
                                                           cmd_path=self.difference_calc_folder_path)
                difference = first_energy - second_energy
                zeroed = difference - self.reference_E
                print('singlet - triplet - reference')
                print('  %s - %s - %s = %s' % (first_energy, second_energy, self.reference_E, zeroed))
                print('singlet - triplet')
                print('  %s - %s = %s' % (first_energy, second_energy, difference))
                return first_energy - second_energy - self.reference_E
            else:
                energy_for_value = self.find_energy_for_value(trial_value,
                                                              orbital_to_optimise,
                                                              energy_type=energy_type,
                                                              cmd_path=self.calc_folder_path)
                zeroed = energy_for_value - self.reference_E
                print('current - reference')
                print('  %s - %s = %s' % (energy_for_value, self.reference_E, zeroed))
                return numpy.abs(energy_for_value - self.reference_E)
        self.optimisation_parameter = optimisation_parameter
        self.reference_E = reference_energy
        print('values:', self.optimisation_parameter, self.reference_E)
        print('Optimising %s over %s, f(%s) = %s and f(%s) = %s' % (
            optimisation_parameter,
            self.trial_range,
            self.trial_range[0],
            find_e_difference(self.trial_range[0]),
            self.trial_range[1],
            find_e_difference(self.trial_range[1]))
        )
        return scipy.optimize.minimize_scalar(find_e_difference)
        # return brent_opt(find_e_difference, self.trial_range[0], self.trial_range[1])
        ## return scipy.optimize.minimize(find_e_difference, numpy.array([self.trial_range[0]]))
        ## return brute(find_e_difference, (self.trial_range[0], self.trial_range[1]))


def optimise_energy():
    optimiser = Optimiser()

    ### Ethene
    s_orbital_to_optimise_s = {'line_type': '$ecp',
                               'basis_ecp_name': 'h ecp-s',
                               'orbital_descriptor': 's-f',
                               'orbital_irrep': '1b1u',
                               'spin': 'mos',
                               'functions_list': [{'coefficient': -7.524, 'r^n': 1, 'exponent': 5.0}]
                               }
    s_orbital_to_optimise_t = {
                               'line_type': '$ecp',
                               'basis_ecp_name': 'h ecp-s',
                               'orbital_descriptor': 's-f',
                               'orbital_irrep': '1b3g',
                               'spin': 'mos',
                               'functions_list': [{'coefficient': -7.524, 'r^n': 1, 'exponent': 5.0}]
                               }
    pz_orbital_to_optimise_s = {
                           'line_type': '$ecp',
                           'basis_ecp_name': 'c ecp-p',
                           'orbital_descriptor': 'p-f',
                           'orbital_irrep': '1b1u',
                           'spin': 'mos',
                           'functions_list': [{'coefficient': -3.88139893705, 'r^n': 1, 'exponent': 0.149275103230995}]
                           }
    pz_orbital_to_optimise_t = {
                           'line_type': '$ecp',
                           'basis_ecp_name': 'c ecp-p',
                           'orbital_descriptor': 'p-f',
                           'orbital_irrep': '1b3g',
                           'spin': 'mos',
                           'functions_list': [{'coefficient': -3.88139893705, 'r^n': 1, 'exponent': 0.149275103230995}]
                           }

    ### Ethane
    c2h6_p_pp_to_optimise = {'line_type': '$ecp',
                             'basis_ecp_name': 'c ecp-p',
                             'orbital_descriptor': 'p-f',
                             'orbital_irrep': '5a',
                             'spin': 'mos',
                             'functions_list': [{'coefficient': -2.8056200103,
                                                 'r^n': 1,
                                                 'exponent': 0.6244618784275526}]
                             }

    c2h6_s_pp_to_optimise = {'line_type': '$ecp',
                             'basis_ecp_name': 'he ecp-s',
                             'orbital_descriptor': 's-f',
                             'orbital_irrep': '5a',
                             'spin': 'mos',
                             'functions_list': [{'coefficient': 0.5,
                                                 'r^n': 1,
                                                 'exponent': 1.5}]
                             }

    # reference_energy = -0.489811  # ethane HOMO
    reference_energy = -13.328  # ethane HOMO
    # optimiser.reference_energy = -0.489811 # ethane HOMO

    basis_information = {
        'basis_name': 'c def2-SV(P)',
        'orbital_to_extract': 'p'
    }

    mos_ref_information = {
        'orbital': '1 b1u',
        'coeff_range': [0, 2]
    }

    optimised_value = optimiser.run_minimise_total_e_difference(
                                                                'coefficient',
                                                                #'exponent',
                                                                orbital_to_optimise=c2h6_s_pp_to_optimise,
                                                                #orbital_to_optimise=c2h6_p_pp_to_optimise,
                                                                energy_type='orbital',
                                                                difference_orbital_to_optimise=None,
                                                                reference_energy=reference_energy,
                                                                basis_info=basis_information,
                                                                mos_info=mos_ref_information,
                                                                use_zeff_calc=False)
    print(optimised_value)
    print('Optimisation converged to energy %s eV, with %s %s' % (reference_energy,
                                                                  optimiser.optimisation_parameter,
                                                                  optimised_value))

optimise_energy()



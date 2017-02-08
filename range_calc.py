import linecache
from chem_lib.optimise import BasisControl
from scipy.optimize import brentq as brent_opt


class Optimiser:
    """Contains functions that let the user optimise pseudo-potential parameters based on reference energy."""
    def __init__(self):
        self.reference_E = None
        self.orbital_to_optimise = None
        self.basis = BasisControl()
        self.basis.variable_file_path = 'basis'
        self.trial_range = (-10.0, 10.0)
        self.spin_indices = {'alpha': -2, 'beta': -1}
        self.optimisation_type = 'coefficient'

    def find_energy_for_value(self, trial_value):
        self.orbital_to_optimise['functions_list'][0][self.optimisation_type] = str(trial_value)
        self.basis.update_variable(self.orbital_to_optimise)
        self.basis.run_ridft(add_to_log=True)

        out_file = open('ridft.log', 'r')
        output_energies = []
        for i, line in enumerate(out_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if split_line[0] == "irrep" and self.orbital_to_optimise['orbital_irrep'] in split_line:
                    energy_line = ' '.join(linecache.getline('ridft.log', i + 3).split()).split()
                    orbital_e = energy_line[split_line.index(self.orbital_to_optimise['orbital_irrep'])]
                    print "new energy %s eV" % orbital_e
                    output_energies.append(orbital_e)
                    linecache.clearcache()

        return float(output_energies[self.spin_indices[self.orbital_to_optimise['spin']]]) - self.reference_E

    def run_optimise(self, optimisation_type, reference_energy):
        self.optimisation_type = optimisation_type
        self.reference_E = reference_energy
        return brent_opt(self.find_energy_for_value, self.trial_range[0], self.trial_range[1])


def optimise_energy():
    optimiser = Optimiser()
    optimiser.orbital_to_optimise = {'line_type': '$ecp',
                                     'basis_ecp_name': 'h ecp-s',
                                     'orbital_descriptor': 's-f',
                                     'orbital_irrep': '1a2\"',
                                     'spin': 'alpha',
                                     'functions_list': [{'coefficient': -7.524, 'r^n': 1, 'exponent': 0.5}]
                                     }
    optimiser.orbital_to_optimise = {'line_type': '$ecp',
                                     'basis_ecp_name': 'c ecp-p',
                                     'orbital_descriptor': 'p-f',
                                     'orbital_irrep': '1a2\"',
                                     'spin': 'alpha',
                                     'functions_list': [{'coefficient': 3.326, 'r^n': 1, 'exponent': 0.295}]
                                     }
    reference_energy = -10.537  # eV
    optimised_value = optimiser.run_optimise('coefficient', reference_energy)
    print 'Optimisation converged to energy %s eV, with %s %s' % (reference_energy,
                                                                  optimiser.optimisation_type,
                                                                  optimised_value)

optimise_energy()


def test_calcs():
    orbital_to_update = {'line_type': '$ecp',
                         'basis_ecp_name': 'h ecp-s',
                         'orbital_descriptor': 's-f'}
    new_functions = [
                    {'exponent': -7.524, 'r^n': 1, 'coefficient': 100.0}
                    ]
    orbital_to_update['functions_list'] = new_functions
    trial_coefficients = [-50.0, -20.0, -10.0, -1.0, -0.1, 0.0, 0.1, 0.2, 0.5, 1.0, 10.0, 100.0, 1000.0, 1200.0, 1500.0]
    trial_coefficients_2 = [250, 500, 750]
    trial_exponents = [100.0, 250.0, 500.0, 550.0, 600.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 750.0, 1000.0, 1500.0, 2000.0]
    trial_coefficients_3 = [x / 1000 for x in range(-3800, -3500, 5)]

    trial_values = trial_exponents

    for coeff in trial_values:
        basis = BasisControl()
        basis.variable_file_path = 'basis (copy)'
        orbital_to_update['functions_list'][0]['coefficient'] = str(coeff)
        basis.update_variable(orbital_to_update)
        basis.run_ridft()

    out_file = open('calcs.txt', 'r')
    output_energies = []
    print 'Checking for energies...'
    for i, line in enumerate(out_file):
        split_line = ' '.join(line.split()).split()
        # if len(split_line) == 6:
        #     if split_line[1] == "total" and split_line[2] == "energy":
        #         total_e = split_line[4]
        #         print trial_values
        #         print "found energy %s" % total_e
        #         output_energies.append(total_e)
        if len(split_line) == 6:
            if split_line[0] == "irrep" and split_line[1] == "1b1u":
                energy_line = ' '.join(linecache.getline('calcs.txt', i+2).split()).split()
                total_e = energy_line[2]
                print trial_values
                print "found energy %s" % total_e
                output_energies.append(total_e)
    print trial_values
    print "energies: %s" % output_energies

    with open('trial_outcomes.txt', 'w') as var_file:
        for index, trial_value in enumerate(trial_values):
            var_file.writelines('%s %s \n' % (trial_value, output_energies[index]))

# test_calcs()

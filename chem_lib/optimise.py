from __future__ import division
from collections import OrderedDict
import subprocess

"""
This library allows the user to read and modify code in a TurboMole basis file.
The user can also reference individual parameters e.g. the coefficient of the
first 2p orbital function in a c TZ basis would be:

coefficient = parsed_basis['$basis']['c TZ']['2 p'][0]['coefficient']

The r^n value of the second 1s orbital function in a h ecp-sp3 ecp would be:

r^n = parsed_basis['$ecp']['h ecp-sp3']['1 s'][1]['r^n']

Example structure of parsed basis file below.
"""

example_basis_file_dict = {
    '$basis': {'h def-SV(P)': {'1 s': [
                                       {'coefficient': 100, 'exponent': 100},
                                       {'coefficient': 200, 'exponent': 200},
                                       ],
                               '2 s': [
                                       {'coefficient': 100, 'exponent': 100}
                                       ]
                               },
               'c def-SV(P)': {'1 s': [
                                       {'coefficient': 100, 'exponent': 100}
                                       ]
                               }
               },
    '$ecp': {'h ecp': {'example (s-f)': [
                                         {'coefficient': 100, 'r^n': 100, 'exponent': 100}
                                         ]
                       }
             }
}


class BasisControl:
    def __init__(self):
        self.variable_file_path = '../tm_files/basis'
        self.line_types = ['$basis', '$ecp']
        self.orbital_descriptors = ['1s', '2s', '3s', '4s', '5s',
                                    '1p', '2p', '3p', '4p', '5p',
                                    '1d', '2d', '3d', '4d', '5d',
                                    'f', 's-f', 'p-f', 'd-f', 'd', 's-d', 'p-d', 'p', 's-p']
        self.basis_file = None

    def read_variables(self):
        variables = {'$basis': OrderedDict(),
                     '$ecp': OrderedDict()}
        line_type = ''
        basis_ecp_name = ''
        orbital_descriptor = ''
        var_file = open(self.variable_file_path, 'r')
        print('Opened file: %s' % self.variable_file_path)

        for line in var_file:
            function_dict = None
            stripped = ''.join(line.split())
            splitted = ' '.join(line.split()).split()

            # Update line and orbital type if needed.
            if stripped in self.line_types:
                line_type = stripped

            elif line[0] != ' ' and len(splitted) is 2:
                basis_ecp_name = ' '.join(line.split())
                if basis_ecp_name not in variables[line_type].keys():
                    print('found basis: %s' % basis_ecp_name)
                    variables[line_type][basis_ecp_name] = OrderedDict()

            elif stripped in self.orbital_descriptors:
                orbital_descriptor = stripped

                if orbital_descriptor not in variables[line_type][basis_ecp_name].keys():
                    variables[line_type][basis_ecp_name][orbital_descriptor] = []
                else:
                    similar_basis_functions = len(["badgers" for key in variables[line_type][basis_ecp_name].keys()
                                                   if orbital_descriptor in key]) + 1
                    variables[line_type][basis_ecp_name][orbital_descriptor + ' (%s)' % similar_basis_functions] = []
                    orbital_descriptor += ' (%s)' % similar_basis_functions

            # Add functions if present in line.
            if len(splitted) in [2, 3]:
                try:
                    if line_type == '$basis':
                        function_dict = {'exponent': float(splitted[0]),
                                         'coefficient': float(splitted[1])}
                    elif line_type == '$ecp':
                        function_dict = {'exponent': float(splitted[0]),
                                         'r^n': int(splitted[1]),
                                         'coefficient': float(splitted[2])}

                except Exception:
                    continue

            if function_dict:
                variables[line_type][basis_ecp_name][orbital_descriptor].append(function_dict)

        var_file.close()
        self.basis_file = variables

    def update_variable(self, new_variable):
        var_file = open(self.variable_file_path, 'r')
        var_file_data = var_file.readlines()

        line_type = ''
        basis_ecp_name = ''
        orbital_descriptor = ''

        # Search for the correct orbital
        for lineno, line in enumerate(var_file_data):
            stripped = "".join(line.split())
            if stripped in self.line_types:
                line_type = stripped
            elif line[0] != ' ' and len((' '.join(line.split())).split()) is 2:
                basis_ecp_name = ' '.join(line.split())
            elif stripped in self.orbital_descriptors:
                orbital_descriptor = stripped

            # Add new details
            if line_type == new_variable['line_type'] \
                    and basis_ecp_name == new_variable['basis_ecp_name'] \
                    and orbital_descriptor == new_variable['orbital_descriptor']:

                for funcno, new_arguments in enumerate(new_variable['functions_list'], start=0):
                    func = new_variable['functions_list'][funcno]

                    if line_type == '$basis':
                        var_file_data[lineno+funcno+1] = ' %s %s \n' % (func['coefficient'],
                                                                        func['exponent'])
                    elif line_type == '$ecp':
                        var_file_data[lineno+funcno+1] = ' %s %s %s \n' % (func['coefficient'],
                                                                           func['r^n'],
                                                                           func['exponent'])

                    # And write everything back
                    with open(self.variable_file_path, 'w') as var_file:
                        if var_file_data:
                            var_file.writelines(var_file_data)

                    print("Modified %s: %s, %s with new coeff %s, exp %s" % (self.variable_file_path,
                                                                             new_variable['basis_ecp_name'],
                                                                             new_variable['orbital_descriptor'],
                                                                             func['coefficient'],
                                                                             func['exponent']))

                var_file.close()
                return

    def get_coefficients_exponents(self, basis_file_extract, basis_set, shell_type):
        """Retrieve coefficients and exponents from basis file."""
        alpha_list = []
        coeff_list = []

        def extract_for_shell(basis_function_name):
            sub_alpha_list = []
            sub_coeff_list = []
            for shell in basis_file_extract["$basis"][basis_set].keys():
                if basis_function_name in shell:
                    a = basis_file_extract["$basis"][basis_set]
                    for exp in a[shell]:
                        sub_alpha_list.append(exp["exponent"])
                        sub_coeff_list.append(exp["coefficient"])
            return sub_alpha_list, sub_coeff_list

        p_alpha, p_coeff = extract_for_shell(shell_type)
        coeff_list.append(p_coeff)
        alpha_list.append(p_alpha)

        coeff_list = [item for sublist in coeff_list for item in sublist]
        alpha_list = [item for sublist in alpha_list for item in sublist]

        return coeff_list, alpha_list

    def extract_contracted_coefficients(self, basis_name, orbital_function_type):
        contracted_coefficients_list = []
        contracted_exponents_list = []
        for function_type, functions in self.basis_file['$basis'][basis_name].items():
            if orbital_function_type in function_type:
                contracted_coefficients_list.append([function['coefficient'] for function in functions])
                contracted_exponents_list.append([function['exponent'] for function in functions])
        return contracted_coefficients_list, contracted_exponents_list

    @staticmethod
    def run_dscf(add_to_log=False, file_path='', silent=False, shell=True, custom_log_path=None):
        if silent is False:
            print("running %s dscf, add to log is %s" % (file_path, add_to_log))
        if custom_log_path:
            command = 'dscf | tee dscf.log %s; cp mos mos_%s' % (custom_log_path, custom_log_path)
        elif add_to_log == True:
            command = 'dscf > dscf.log'
        else:
            command = 'dscf'
        subprocess.call(command, shell=shell, cwd=file_path)

    @staticmethod
    def run_escf(add_to_log=False, file_path='', silent=False, shell=True):
        if silent is False:
            print("running %s escf, add to log is %s" % (file_path, add_to_log))
        if add_to_log == True:
            command = 'escf > escf.log'
        else:
            command = 'escf'
        subprocess.call(command, shell=shell, cwd=file_path)

    @staticmethod
    def run_ridft(add_to_log=False, file_path='', silent=False, shell=True):
        if silent is False:
            print("running %s ridft, add to log is %s" % (file_path, add_to_log))
        if add_to_log == True:
            command = 'ridft > ridft.log'
        else:
            command = 'ridft'
        subprocess.call(command, shell=shell)

    @staticmethod
    def run_eiger(add_to_log=False, file_path='', silent=False, shell=True):
        if silent is False:
            print("running %s eiger, add to log is %s" % (file_path, add_to_log))
        if add_to_log == True:
            command = 'eiger > eiger.log'
        else:
            command = 'eiger'
        subprocess.call(command, shell=shell)

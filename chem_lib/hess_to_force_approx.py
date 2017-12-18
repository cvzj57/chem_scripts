"""
This script converts hessian format and force format turbomole files.
"""


class HessianControl:
    def __init__(self):
        self.variable_file_path = '../tm_files/'
        self.max_elements_per_line = 8

    def generate_forceapprox_file_from_hessian(self):
        """Reads Hessian into list."""
        var_file = open(self.variable_file_path+'hessapprox', 'r')
        matrix_element_list = []

        for line in var_file:
            splitted = line.split()

            if splitted[0] != '$hessapprox' and splitted[0] != '$end':
                matrix_element_list += splitted

        print('Diagonal matrix size: ' + str(len(matrix_element_list)))

        target_elements_this_line = 1
        actual_elements_this_line = 0
        new_force_matrix = []
        current_line = []

        for matrix_element in matrix_element_list:
            print(current_line)
            if actual_elements_this_line == target_elements_this_line:
                target_elements_this_line += 1
                actual_elements_this_line = 1
                new_force_matrix.append(current_line)
                current_line = [float(matrix_element.replace('D', 'e'))]
            # elif len(current_line) == self.max_elements_per_line:
            #     actual_elements_this_line += 1
            #     new_force_matrix.append(current_line)
            #     current_line = [float(matrix_element.replace('D', 'e'))]
            else:
                # current_line.insert(-1, float(matrix_element.replace('D', 'e')))
                current_line.append(float(matrix_element.replace('D', 'e')))
                actual_elements_this_line += 1
            print(current_line)

        print('Length of forceapprox matrix: ' + str(len([item for sublist in new_force_matrix for item in sublist])))

        new_force_file = open(self.variable_file_path+'forceapprox', 'w+')
        new_force_file.writelines('$forceapprox\n')
        for line in new_force_matrix:
            new_force_file.writelines(' '+' '.join([str((" % .5f" % element)) for element in line])+'\n')
        new_force_file.writelines('$end\n')
        new_force_file.close()

        print('Generated forceapprox file from hessapprox matrix.')


if __name__ == "__main__":
    control = HessianControl()
    control.generate_forceapprox_file_from_hessian()
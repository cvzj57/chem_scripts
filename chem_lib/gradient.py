"""
This library contains code for handling Turbomole's gradient files.
"""

sample_cycle = [
    {'cycle': 1,
     'energy': -40.0,
     'dE': -0.05,
     'atoms': [{'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 1,
               'dx': -0.5, 'dy': 0.5, 'dz': 0.0}]
     },
]


class GradientControl:
    def __init__(self):
        self.variable_file_path = '../tm_files/gradient'
        self.gradient_file = []

    def read_variables(self):
        """Reads Turbomole gradients into init variables."""
        var_file = open(self.variable_file_path, 'r')
        print('Opened file: %s' % self.variable_file_path)

        cycle_dict = {'atoms': []}
        current_coord_hash = 1
        current_gradient_hash = 1

        for line in var_file:
            splitted = line.split()

            if splitted[0] == 'cycle':
                cycle_dict['cycle'] = splitted[2]
                cycle_dict['energy'] = splitted[6]
                cycle_dict['dE'] = splitted[9]
                current_coord_hash = 1
                current_gradient_hash = 1
                if len(cycle_dict['atoms']) != 0:
                    self.gradient_file.append(cycle_dict)

            if len(splitted) == 4 and (len(splitted[-1]) == 1 or len(splitted[-1]) ==2):
                cycle_dict['atoms'].append({'x': splitted[0],
                                            'y': splitted[1],
                                            'z': splitted[2],
                                            'el': splitted[-1],
                                            '#': current_coord_hash})
                current_coord_hash += 1

            if len(splitted) == 3 and line[0] != '$':
                atom_index = next(index for (index, d) in enumerate(cycle_dict['atoms']) if d['#'] == current_gradient_hash)
                cycle_dict['atoms'][atom_index]['dx'] = float(splitted[0].replace('D', 'E'))
                cycle_dict['atoms'][atom_index]['dy'] = float(splitted[1].replace('D', 'E'))
                cycle_dict['atoms'][atom_index]['dz'] = float(splitted[2].replace('D', 'E'))
                current_gradient_hash += 1

        self.gradient_file.append(cycle_dict)

        print('Successfully read gradient file (%s cycle(s), %s atom(s))' % (len(self.gradient_file),
                                                                             len(self.gradient_file[0]['atoms'])
                                                                             ))

    def translate_into_fortran(self, dEdZ):
        new_dE_string = ("%e" % dEdZ).replace('e', 'D')
        if dEdZ[0] == '-' and dEdZ[2] == '.':
            if dEdZ[-1] == '0':
                new_ending = 1
            elif dEdZ[-3] == '+':
                new_ending = int(dEdZ[-1]) + 1
            elif dEdZ[-3] == '-':
                new_ending = int(dEdZ[-1]) - 1
            new_beginning = dEdZ[0] + '.' + dEdZ[1]
            new_dE_string = new_beginning + dEdZ[3:-1] + str(new_ending)
        return new_dE_string


if __name__ == "__main__":
    control = GradientControl()
    control.read_variables()
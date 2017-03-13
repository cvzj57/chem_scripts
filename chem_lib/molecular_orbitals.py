from collections import OrderedDict

"""
This library contains code for handling Turbomole's molecular orbital files.
"""

example_mos_dict = {
    '1 a1': ['0.1234567890123D+12', '0.34567890123D-08'],
    '4 b2': ['-.9876543219876D+14', '0.76543219876D+02']
    }


class MOSControl:
    def __init__(self):
        self.variable_file_path = 'alpha'
        self.file_description = None
        self.values_per_line = None
        self.chars_per_value = None
        self.molecular_orbital_file = OrderedDict()

    def read_variables(self):
        """Reads Turbomole orbitals into init variables."""
        var_file = open(self.variable_file_path, 'r')
        print('Opened file: %s' % self.variable_file_path)

        current_orbital = None
        current_values = []

        for line in var_file:
            stripped = ''.join(line.split())
            splitted = (' '.join(line.split())).split()

            if stripped[0] == '$' and stripped != '$end':
                print('Found file description %s' % ' '.join(line.split()))
                self.file_description = stripped
                line_format = stripped[stripped.find('format(')+7:stripped.index(')')]
                self.values_per_line = int(line_format.split('d')[0])
                self.chars_per_value = int(line_format.split('d')[1].split('.')[0])

            if len(splitted) == 4:
                if current_orbital:
                    self.molecular_orbital_file[current_orbital] = [float(value.replace('D', 'E').replace('-.', '-0.'))
                                                                    for sublist in current_values
                                                                    for value in sublist if value]
                current_orbital = splitted[0]+' '+splitted[1]
                current_values = []

            if stripped[0] == '-' or stripped[0] == '0':
                current_values.append([stripped[i:i+self.chars_per_value]
                                       for i in range(0,
                                                      self.values_per_line*self.chars_per_value,
                                                      self.chars_per_value)])

        self.molecular_orbital_file[current_orbital] = [float(value.replace('D', 'E').replace('-.', '-0.'))
                                                        for sublist in current_values for value in sublist if value]

        print('Found MOs: %s' % ', '.join([key for key in self.molecular_orbital_file.keys()]))

if __name__ == "__main__":
    control = MOSControl()
    control.read_variables()

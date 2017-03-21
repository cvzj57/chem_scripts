import numpy
import math

sample_coords = [
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 1},
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'h', '#': 2}
]


class CoordControl:
    def __init__(self):
        self.coord_file_path = '../tm_files/coord_c4h6'
        self.coord_list = []

    def read_coords(self):
        var_file = open(self.coord_file_path, 'r')

        for lineno, line in enumerate(var_file):
            stripped = ''.join(line.split())
            splitted = ' '.join(line.split()).split()

            new_atom = {}

            if stripped[0] == '$':
                if 'coord' in stripped:
                    print('Found coords.')
                    continue
                else:
                    print('No more coords.')
                    break

            new_atom['el'] = splitted[-1]
            new_atom['z'] = float(splitted[-2])
            new_atom['y'] = float(splitted[-3])
            new_atom['x'] = float(splitted[-4])
            new_atom['#'] = lineno

            self.coord_list.append(new_atom)
        print(self.coord_list)

    def write_coord(self, new_atom):
        var_file = open(self.coord_file_path, 'r')
        var_file_data = var_file.readlines()

        var_file_data[new_atom['#']] = ' %s %s %s %s\n' % (
            new_atom['x'],
            new_atom['y'],
            new_atom['z'],
            new_atom['el'])

        # And write everything back
        with open(self.coord_file_path, 'w') as var_file:
            if var_file_data:
                var_file.writelines(var_file_data)

        print("Modified atom %s, now %s %s %s %s" % (
            new_atom['#'],
            new_atom['x'],
            new_atom['y'],
            new_atom['z'],
            new_atom['el']))

        var_file.close()

    def delete_hydrogen_atoms(self):
        var_file = open(self.coord_file_path, 'r')
        var_file_data = var_file.readlines()

        for line in var_file_data:
            print('|'+line+'|')
        #
        # for lineno, line in enumerate(var_file_data):
        #     stripped = ''.join(line.split())
        #     splitted = ' '.join(line.split()).split()
        #
        #     print(lineno)
        #     print(stripped)
        #
        #     if line[0] == '$':
        #         if '$coord' in stripped:
        #             print('Found coords.')
        #             continue
        #         else:
        #             print('No more coords.')
        #             break
        #
        #     if 'h' in line:
        #         del var_file_data[lineno]
        #
        # # And write everything back
        # with open(self.coord_file_path, 'w') as var_file:
        #     if var_file_data:
        #         var_file.writelines(var_file_data)
        #
        # print("Removed all hydrogen atoms.")
        #
        # var_file.close()

    def vectorise_atom(self, index):
        return numpy.array([float(self.coord_list[index-1]['x']),
                            float(self.coord_list[index-1]['y']),
                            float(self.coord_list[index-1]['z'])])

    def measure_atom_atom_dist(self, index_1, index_2):
        return numpy.linalg.norm(self.vectorise_atom(index_2) - self.vectorise_atom(index_1))

    def atom_atom_vector(self, index_1, index_2):
        return self.vectorise_atom(index_2) - self.vectorise_atom(index_1)

    def find_nearest_carbon(self, atom_index):

        carbon_indices = []
        distances = []
        for atom in self.coord_list:
            if atom_index != atom['#']:
                if atom['el'] == 'c':
                    carbon_indices.append(atom['#'])
                    distances.append(control.measure_atom_atom_dist(atom_index, atom['#']))

        index = carbon_indices[numpy.argmin(distances)]
        return index, numpy.min(distances)

    def order_atoms_by_distance_from(self, central_atom_index):

        def distance_from(list_atom):
            return self.measure_atom_atom_dist(central_atom_index, list_atom['#'])

        return sorted(self.coord_list, key=distance_from)

    def lengtherise_vector(self, vector, target_length):
        current_length = numpy.linalg.norm(vector)
        print(current_length)
        return vector * (target_length/current_length)

    def construct_euler_rodriguez_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = numpy.asarray(axis)
        axis = axis / math.sqrt(numpy.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return numpy.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                           [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                           [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


if __name__ == "__main__":
    control = CoordControl()
    control.read_coords()
    # control.write_coord({'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 1})
    practice_vector = numpy.array([2.0, 10.0, 1.0])
    print(control.lengtherise_vector(practice_vector, 0.5))
    print(numpy.linalg.norm(control.lengtherise_vector(practice_vector, 0.5)))
    # control.delete_hydrogen_atoms()

    vector_1 = control.atom_atom_vector(1, 2)
    vector_2 = control.atom_atom_vector(1, 3)

    print(vector_1, vector_2)

    distanced_atom_list = control.order_atoms_by_distance_from(1)

    print(numpy.cross(control.vectorise_atom(distanced_atom_list[1]['#']), control.vectorise_atom(distanced_atom_list[2]['#'])))


import numpy
import math
import time

sample_coords = [
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 1},
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'h', '#': 2}
]


class CoordControl:
    def __init__(self):
        self.coord_file_path = '../tm_files/coord_c4h6'
        self.coord_list = []
        self.no_potential_sets_per_atom = 3
        self.atom_potential_set_distance = 0.5
        self.potential_set_split_distance = 0.25

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
        var_file.close()

    def write_coord(self, new_atom, overwrite=False):
        var_file = open(self.coord_file_path, 'r')
        var_file_data = var_file.readlines()
        self.coord_list = []
        self.read_coords()

        def return_hash(atom): return atom['#']

        if overwrite is True:
            atom_hash = new_atom['#']
            del var_file_data[atom_hash]
        else:
            print(sorted(self.coord_list, key=return_hash, reverse=True)[0]['#'])
            atom_hash = sorted(self.coord_list, key=return_hash, reverse=True)[0]['#']+1

        var_file_data.insert(atom_hash, ' %s %s %s %s\n' % (
            new_atom['x'],
            new_atom['y'],
            new_atom['z'],
            new_atom['el'])
        )

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
        lines_to_delete = []

        for lineno, line in enumerate(var_file_data):
            stripped = ''.join(line.split())

            if line[0] == '$':
                if '$coord' in stripped:
                    print('Found coords.')
                    continue
                else:
                    print('No more coords.')
                    break

            if 'h' in line:
                lines_to_delete.insert(0, lineno)

        for lineno in lines_to_delete:
            del var_file_data[lineno]

        # And write everything back
        with open(self.coord_file_path, 'w') as var_file:
            if var_file_data:
                var_file.writelines(var_file_data)

        print("Removed all hydrogen atoms.")

        var_file.close()

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

    def order_atoms_by_distance_from(self, central_atom_index, element=None):

        if element:
            coord_list = (item for item in self.coord_list if item["el"] == element)
        else:
            coord_list = self.coord_list

        def distance_from(list_atom):
            # print(self.measure_atom_atom_dist(central_atom_index, list_atom['#']))
            return self.measure_atom_atom_dist(central_atom_index, list_atom['#'])

        return sorted(coord_list, key=distance_from)

    def lengtherise_vector(self, vector, target_length):
        current_length = numpy.linalg.norm(vector)
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

    def pseudopotentialise_molecule(self):

        atoms_to_potentialise = (item for item in self.coord_list if item["el"] == 'c')
        potential_coords_list = []

        for atom in atoms_to_potentialise:
            distanced_atom_list = self.order_atoms_by_distance_from(atom['#'])
            distanced_carbon_list = self.order_atoms_by_distance_from(atom['#'], element='c')
            # for carbon in distanced_carbon_list:
            #     print(carbon)
            print(distanced_carbon_list[1]['#'], atom['#'])
            primary_vector = self.vectorise_atom(distanced_carbon_list[1]['#'])-self.vectorise_atom(atom['#'])
            print(numpy.linalg.norm(primary_vector))
            normal_vector = numpy.cross(
                self.vectorise_atom(distanced_atom_list[1]['#']) - self.vectorise_atom(atom['#']),
                self.vectorise_atom(distanced_atom_list[2]['#']) - self.vectorise_atom(atom['#'])
            )

            primary_potential_vector = self.lengtherise_vector(primary_vector, self.atom_potential_set_distance)
            # print(numpy.linalg.norm(primary_potential_vector))
            potential_set_split_vector = self.lengtherise_vector(normal_vector, self.potential_set_split_distance)

            relative_potential_vectors = [
                primary_potential_vector + potential_set_split_vector,
                primary_potential_vector - potential_set_split_vector
            ]

            for potential_set in range(self.no_potential_sets_per_atom-1):

                pps_positive = numpy.dot(self.construct_euler_rodriguez_matrix(
                                             normal_vector,
                                             2*numpy.pi/self.no_potential_sets_per_atom),
                                         relative_potential_vectors[-2],
                                         )
                pps_negative = numpy.dot(self.construct_euler_rodriguez_matrix(
                                             normal_vector,
                                             2*numpy.pi/self.no_potential_sets_per_atom),
                                         relative_potential_vectors[-1]
                                         )

                relative_potential_vectors.append(pps_positive)
                relative_potential_vectors.append(pps_negative)


            # potential coords are still relative to their atom, now make them real
            for vector in relative_potential_vectors:
                # print({'#': atom['#'], 'el': 'h', 'x': vector[0]+atom['x'], 'y': vector[0]+atom['y'], 'z': vector[0]+atom['z']})
                potential_coords_list.append(
                    {'#': 0, 'el': 'h', 'x': vector[0]+atom['x'], 'y': vector[1]+atom['y'], 'z': vector[2]+atom['z']},
                )

        # print(potential_coords_list)

        # Now add potentials to coord list, after removing the 'real' hydrogen atoms.
        self.delete_hydrogen_atoms()
        time.sleep(5)
        print("here come dem potentials")
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=False)




if __name__ == "__main__":
    control = CoordControl()
    control.read_coords()
    # control.write_coord({'x': 1.24869729097288, 'y': -0.35894369997066, 'z': -0.48876369930179, 'el': 'c', '#': 1}, overwrite=True)
    # control.write_coord({'x': -1.24781026329849, 'y': 0.35940201268517, 'z': 0.49082623491110, 'el': 'c', '#': 6})
    # control.write_coord({'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 3})

    print("potentialising")
    control.pseudopotentialise_molecule()


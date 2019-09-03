import numpy
import math
import sys
import os
import subprocess
import shutil
'''Library for the reading and manipulation of turbomole coordinate files. Contains functions to supply
pseudo-potentials to molecules.'''

sample_coords = [
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'c', '#': 1},
    {'x': 0.0, 'y': 0.0, 'z': 0.0, 'el': 'h', '#': 2}
]

define_cmds = '''

a coord
*
no
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
*



q


q
'''


class CoordControl:
    def __init__(self):
        self.coord_file_path = './coord'
        self.coord_list = []
        self.no_potential_sets_per_atom = 3
        self.sp2_1e_atom_potential_set_distance = 0.58206525
        self.sp2_2e_atom_potential_set_distance = 0.5
        self.sp3_1e_atom_potential_distance = 0.5
        self.sp2_1e_potential_set_split_distance = 0.26891659
        self.sp2_2e_potential_set_split_distance = 0.25
        self.sp3_pseudo_element = 'li'
        self.sp3_2e_pseudo_element = 'li'
        self.sp2_pseudo_element = 'he'
        self.all_e_basis = 'def-SV(P)'
        self.sp2_pseudo_carbon_basis = 'sless-SV(P)'
        self.sp3_pseudo_carbon_basis = 'def-SV(P)'
        self.sp3_2e_pseudo_carbon_basis = 'def-SV(P)'
        self.sp3_carbon_ecp = 'sp3-ecp-p'
        self.sp3_2e_carbon_ecp = 'sp3-ecp-p2'
        self.sp2_carbon_ecp = 'sp2-ecp-p'
        self.sp2_hydrogen_ecp = 'sp2-ecp-s'
        self.sp3_hydrogen_ecp = 'sp3-ecp-s'
        self.sp3_2e_hydrogen_ecp = 'sp3-ecp-s2'
        self.sp2_2e_carbon_ecp = 'sp2-ecp-p2'
        self.sp2_2e_hydrogen_ecp = 'sp2-ecp-s2'
        self.add_primary_vector_potentials_as_coords = True
        self.bond_deciding_distance = 3.7
        self.pseudo_deciding_distance = 1.5
        self.reference_guess_basis_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../tm_files/reference_guess_basis')

    @staticmethod
    def check_is_int(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    def run_define(self, cmds):
        command = 'define < %s' % cmds
        print('Running define...')
        subprocess.call(command, shell=True, cwd=os.path.dirname(os.path.abspath(self.coord_file_path)))

    def read_coords(self):
        var_file = open(self.coord_file_path, 'r')

        for lineno, line in enumerate(var_file):
            stripped = ''.join(line.split())
            splitted = ' '.join(line.split()).split()

            new_atom = {}

            if stripped[0] == '$':
                if 'coord' in stripped:
                    continue
                else:
                    break

            new_atom['el'] = splitted[-1]
            new_atom['z'] = float(splitted[-2])
            new_atom['y'] = float(splitted[-3])
            new_atom['x'] = float(splitted[-4])
            new_atom['#'] = lineno

            self.coord_list.append(new_atom)
        # print('Read atom coords.')
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

        # print("Modified atom %s, now %s %s %s %s" % (
        #     atom_hash,
        #     new_atom['x'],
        #     new_atom['y'],
        #     new_atom['z'],
        #     new_atom['el']))

        var_file.close()

    def parse_coord_list(self, raw_input):
        print("Raw input received: %s" % raw_input)
        splitted = raw_input.split(',')
        final_coord_list = []
        for split_input in splitted:
            if '-' in split_input:
                input_range = split_input.split('-')
                final_coord_list.extend([*range(int(input_range[0]), int(input_range[1])+1)])
            else:
                final_coord_list.append(int(split_input))
        #print(final_coord_list)
        return final_coord_list

    def delete_hydrogen_atoms(self):
        var_file = open(self.coord_file_path, 'r')
        var_file_data = var_file.readlines()
        lines_to_delete = []

        for lineno, line in enumerate(var_file_data):
            stripped = ''.join(line.split())

            if line[0] == '$':
                if '$coord' in stripped:
                    continue
                else:
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

    def delete_specified_atoms(self, deletion_hash_list):
        var_file = open(self.coord_file_path, 'r')
        var_file_data = var_file.readlines()
        lines_to_delete = []

        for lineno, line in enumerate(var_file_data):
            stripped = ''.join(line.split())

            if line[0] == '$':
                if '$coord' in stripped:
                    continue
                else:
                    break

            if lineno in deletion_hash_list:
                lines_to_delete.insert(0, lineno)

        for lineno in lines_to_delete:
            del var_file_data[lineno]

        # And write everything back
        with open(self.coord_file_path, 'w') as var_file:
            if var_file_data:
                var_file.writelines(var_file_data)

        print("Removed atoms: %s" % deletion_hash_list)
        var_file.close()

    def vectorise_atom(self, index):
        return numpy.array([float(self.coord_list[index-1]['x']),
                            float(self.coord_list[index-1]['y']),
                            float(self.coord_list[index-1]['z'])])

    def measure_atom_atom_dist(self, index_1, index_2):
        return numpy.linalg.norm(self.vectorise_atom(index_2) - self.vectorise_atom(index_1))

    def distance_from(self, this_atom, list_atom):
        return self.measure_atom_atom_dist(this_atom['#'], list_atom['#'])

    def order_atoms_by_distance_from(self, central_atom_index, element=None, list_of_atoms_to_distance=None):

        if element == '!h':
            coord_list = (item for item in self.coord_list if item["el"] != 'h')
        elif element == '!pseudo':
            coord_list = (item for item in self.coord_list if item['el'] != self.sp2_pseudo_element and item['el'] != self.sp3_pseudo_element)
        elif element == '!pseudo!h':
            coord_list = (item for item in self.coord_list if item['el'] != self.sp2_pseudo_element and item['el'] != self.sp3_pseudo_element and item['el'] != 'h')
        elif element:
            coord_list = (item for item in self.coord_list if item["el"] == element)
        elif list_of_atoms_to_distance:
            coord_list = list_of_atoms_to_distance
        else:
            coord_list = self.coord_list

        def distance_from(list_atom):
            return self.measure_atom_atom_dist(central_atom_index, list_atom['#'])

        return sorted(coord_list, key=distance_from)

    def lengtherise_vector(self, vector, target_length):
        current_length = numpy.linalg.norm(vector)
        return vector * (float(target_length)/float(current_length))

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

    def get_atoms_within_bond_distance(self, atom_hash):
        distanced_list = self.order_atoms_by_distance_from(atom_hash)
        for distanced_atom in distanced_list:
            if self.measure_atom_atom_dist(atom_hash, distanced_atom['#']) < self.bond_deciding_distance:
                return distanced_list[0:distanced_list.index(distanced_atom)]

    def set_potential_distance_to(self, atom_hash, new_distance):
        """Moves potentials around an atom to specified distance."""
        # identify potentials from atom hash
        pseudopotentials = self.identify_pseudocarbon_potentials(atom_hash)
        potential_coords_list = []
        deletion_list = []

        # Get vectors from potentials to their atom
        # extend to new_distance
        for pseudopotential in pseudopotentials:
            vector_from_pseudo_carbon = self.vectorise_atom(pseudopotential['#']) - self.vectorise_atom(atom_hash)
            new_vector_from_pseudocarbon = self.lengtherise_vector(vector_from_pseudo_carbon, new_distance)
            new_potential_coordinates = self.vectorise_atom(atom_hash) + new_vector_from_pseudocarbon
            potential_coords_list.append(
                    {'#': pseudopotential['#'],
                     'el': self.sp3_pseudo_element,
                     'x': new_potential_coordinates[0],
                     'y': new_potential_coordinates[1],
                     'z': new_potential_coordinates[2]},
                )
            deletion_list.append(pseudopotential['#'])

        # overwrite old potentials
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=True)

    def repseudopotentialise_sp2_atom(self, atom_hash, new_set_distance, new_set_split_distance, supplied_normal_vector=None):
        # identify atom from hash
        atom_to_repotentialise = self.vectorise_atom(atom_hash)
        # identify current pseudopotentials
        current_potentials = self.identify_pseudocarbon_potentials(atom_hash)
        def return_hash(atom): return atom['#']
        hashed_potential_list = sorted(current_potentials, key=return_hash)

        # identify a primary vector
        #distanced_carbon_list = self.order_atoms_by_distance_from(atom_hash, element='c')
        distanced_nonpsnonh_list = self.order_atoms_by_distance_from(atom_hash, element='!pseudo!h')
        #primary_vector = self.vectorise_atom(distanced_carbon_list[1]['#']) - self.vectorise_atom(atom_hash)
        print('Primary vector is to a %s atom. Hope that\'s right.' % distanced_nonpsnonh_list[1]['el'])
        primary_vector = self.vectorise_atom(distanced_nonpsnonh_list[1]['#']) - self.vectorise_atom(atom_hash)

        # new potentials
        new_potential_coords_list = []

        # Identify normal vector to sp2 plane, take any potential, find its partner, and get vector between them
        distanced_potential_list = self.order_atoms_by_distance_from(current_potentials[0]['#'],
                                                                     list_of_atoms_to_distance=current_potentials)

        if supplied_normal_vector is not None:
            normal_vector = supplied_normal_vector
        else:
            normal_vector = self.vectorise_atom(hashed_potential_list[0]['#']) - self.vectorise_atom(hashed_potential_list[1]['#'])

        primary_potential_vector = self.lengtherise_vector(primary_vector, float(new_set_distance))
        potential_set_split_vector = self.lengtherise_vector(normal_vector, float(new_set_split_distance))

        # Keep track of relative vectors for working out potentials
        relative_potential_vectors = [
            primary_potential_vector + potential_set_split_vector,
            primary_potential_vector - potential_set_split_vector
        ]

        # Final ones to be added (won't include the primary relative vectors if beta potential)
        final_relative_potential_vectors = []

        # Is it an alpha or beta potential?
        if len(current_potentials) == 4:
            pass
        else:
            final_relative_potential_vectors = [
                primary_potential_vector + potential_set_split_vector,
                primary_potential_vector - potential_set_split_vector
            ]

        # Build up a list of relative potential vectors of the form [primary+ve, primary-ve, secondary+ve, secondary-ve ...]
        for potential_set in range(self.no_potential_sets_per_atom - 1):
            pps_positive = numpy.dot(self.construct_euler_rodriguez_matrix(
                normal_vector,
                2 * numpy.pi / self.no_potential_sets_per_atom),
                relative_potential_vectors[-2],
            )
            pps_negative = numpy.dot(self.construct_euler_rodriguez_matrix(
                normal_vector,
                2 * numpy.pi / self.no_potential_sets_per_atom),
                relative_potential_vectors[-1]
            )

            relative_potential_vectors.append(pps_positive)
            relative_potential_vectors.append(pps_negative)

            final_relative_potential_vectors.append(pps_positive)
            final_relative_potential_vectors.append(pps_negative)

        # potential coords are still relative to their atom, now make them real.
        for i, vector in enumerate(final_relative_potential_vectors):
            new_potential_coords_list.append(
                {'#': current_potentials[i]['#'],
                 'el': self.sp2_pseudo_element,
                 'x': vector[0] + atom_to_repotentialise[0],
                 'y': vector[1] + atom_to_repotentialise[1],
                 'z': vector[2] + atom_to_repotentialise[2]},
            )

        # Now add potentials to coord list, we must overwrite the original potentials using their hashes.
        # This stops the ECP and basis assignments being invalidated.
        for potential_coord in new_potential_coords_list:
            self.write_coord(potential_coord, overwrite=True)

        #print('Re-potentialised atom %s, with set distance %s, split %s' % (atom_hash,
         #                                                                   new_set_distance,
          #                                                                  new_set_split_distance))

    def set_potential_aperture_angle_to(self, atom_hash, new_distance):
        """Moves potentials around an atom to specified distance."""
        #TODO: Finish this.
        pseudopotentials = self.identify_pseudocarbon_potentials(atom_hash)
        potential_coords_list = []
        deletion_list = []

        for pseudopotential in pseudopotentials:
            # get rotation axis via cross-products
            # if 3 atoms within pseudo-distance this is an sp3 pseudo-carbon
            if len(pseudopotentials) == 3:
                pass

            # if 4 atoms within pseudo-distance this is an sp2 2e pseudo-carbon
            elif len(pseudopotentials) == 4:
                pass


            # if 6 atoms within pseudo-distance this is an sp2 pseudo-carbon
            elif len(pseudopotentials) == 6:
                pass

            # apply euler-rodriguez

            vector_from_pseudo_carbon = self.vectorise_atom(pseudopotential['#']) - self.vectorise_atom(atom_hash)
            new_vector_from_pseudocarbon = self.lengtherise_vector(vector_from_pseudo_carbon, new_distance)
            new_potential_coordinates = self.vectorise_atom(atom_hash) + new_vector_from_pseudocarbon

            potential_coords_list.append(new_potential_coordinates)
            deletion_list.append(pseudopotential['#'])

        self.delete_specified_atoms(deletion_list)
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=False)

    def identify_pseudocarbon_potentials(self, atom_hash):
        pseudos = []

        distanced_atoms = self.order_atoms_by_distance_from(atom_hash)
        nearest_6_distances = [self.measure_atom_atom_dist(atom_hash, distanced_atom['#']) for distanced_atom in
                               distanced_atoms[1:7]]
        pseudo_distances = [less_than_distance for less_than_distance in nearest_6_distances if
                            less_than_distance < self.pseudo_deciding_distance]

        # if 3 atoms within pseudo-distance this is an sp3 pseudo-carbon
        if len(pseudo_distances) == 3:
            pseudos = distanced_atoms[1:4]

        # if 2 atoms within pseudo-distance this is an sp3 2e pseudo-carbon
        if len(pseudo_distances) == 2:
            pseudos = distanced_atoms[1:3]

        # if 4 atoms within pseudo-distance this is an sp2 2e pseudo-carbon
        elif len(pseudo_distances) == 4:
            pseudos = distanced_atoms[1:5]

        # if 6 atoms within pseudo-distance this is an sp2 pseudo-carbon
        elif len(pseudo_distances) == 6:
            pseudos = distanced_atoms[1:7]

        print("Identified %s as a %s pseudoed carbon." % (str(atom_hash), str(len(pseudo_distances))))
        return pseudos

    def chunks(self, list_to_chunk, size):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(list_to_chunk), size):
            yield list_to_chunk[i:i + size]

    def supply_ecps_bases_to_define(self, atom_specifier, ecp_basis_type, ecp_basis_name):
        '''Construct basis and ecp assignment blocks for define.'''
        string_to_return = ''

        if type(atom_specifier) is list:
            chunked_indices = self.chunks(atom_specifier, 20)
            for chunk in chunked_indices:
                string_to_return += """%s\n%s\n%s\nfile\nbasis\n""" % \
                        (ecp_basis_type, ','.join(str(index) for index in chunk), ecp_basis_name)
        elif atom_specifier == 'all':
            string_to_return += """%s\n"%s"\n%s\n\n""" % \
                        (ecp_basis_type, atom_specifier, ecp_basis_name)
        elif type(atom_specifier) is str:
            string_to_return += """%s\n"%s"\n%s\nfile\nbasis\n""" % \
                        (ecp_basis_type, atom_specifier, ecp_basis_name)
        return string_to_return

    def guess_potentialisation(self, sysargs):
        """Guesses molecular potentialisation for sp2 and sp3 carbons."""

        print("Guessing potentialisation...")
        print("Copying reference basis...")
        shutil.copyfile(self.reference_guess_basis_path, os.path.join(os.path.dirname(self.coord_file_path), 'basis'))

        sp2_1e_replacement_list = []
        sp2_2e_replacement_list = []
        sp2_deletion_list = []
        sp3_replacement_list = []
        sp3_deletion_list = []
        sp3_2e_replacement_list = []
        sp3_2e_deletion_list = []
        # 'incl' lets us guess the potentialisation of supplied indices
        # 'excl' will guess potentialisation of all carbon atoms except those specified
        if 'incl' in sysargs:
            include_list = self.parse_coord_list(sysargs[sysargs.index('incl')+1])
            carbon_atoms = [atom for atom in self.coord_list if atom["#"] in include_list]
        else:
            carbon_atoms = [atom for atom in self.coord_list if atom["el"] == 'c']
        if 'excl' in sysargs:
            exclude_list = self.parse_coord_list(sysargs[sysargs.index('excl')+1])
            carbon_atoms = [carbon_atom for carbon_atom in carbon_atoms if carbon_atom['#'] not in exclude_list]

        # Sort through carbons to decide what needs potentialising. Find atoms bonded to each carbon
        for atom in carbon_atoms:
            distanced_atoms = self.order_atoms_by_distance_from(atom['#'])
            nearest_4_distances = [self.measure_atom_atom_dist(atom['#'], distanced_atom['#']) for distanced_atom in
                                   distanced_atoms[1:5]]
            bonded_distances = [less_than_distance for less_than_distance in nearest_4_distances if
                                less_than_distance < self.bond_deciding_distance]

            # if 3 bonded atoms, may be sp2, check if they're hydrogens
            if len(bonded_distances) == 3:
                hydrogens_bonded_to_this_atom = [distanced_atom for distanced_atom in distanced_atoms[1:5] if
                                                 distanced_atom['el'] == 'h' and self.measure_atom_atom_dist(atom['#'], distanced_atom['#']) < self.bond_deciding_distance]
                sp2_deletion_list.extend([hydrogen['#'] for hydrogen in hydrogens_bonded_to_this_atom])
                sp2_1e_replacement_list.append(str(atom['#']))

            # if 4 bonded atoms, may be sp3, check if they're hydrogens
            elif len(bonded_distances) == 4:
                hydrogens_bonded_to_this_atom = [distanced_atom for distanced_atom in distanced_atoms[1:5] if
                                                 distanced_atom['el'] == 'h' and self.measure_atom_atom_dist(atom['#'], distanced_atom['#']) < self.bond_deciding_distance]
                if len(hydrogens_bonded_to_this_atom) == 3:
                    sp3_replacement_list.extend([str(hydrogen['#']) for hydrogen in hydrogens_bonded_to_this_atom])
                    sp3_deletion_list.extend([hydrogen['#'] for hydrogen in hydrogens_bonded_to_this_atom])
                elif len(hydrogens_bonded_to_this_atom) == 2:
                    sp3_2e_replacement_list.extend([str(hydrogen['#']) for hydrogen in hydrogens_bonded_to_this_atom])
                    sp3_2e_deletion_list.extend([hydrogen['#'] for hydrogen in hydrogens_bonded_to_this_atom])

        # Now we need to identify the carbons that are going to be 2e sp2 pseudocarbons
        for carbon_hash in sp2_1e_replacement_list:
            distanced_nonhydrogen_list = self.order_atoms_by_distance_from(int(carbon_hash), element='!h')
            nonhydrogens_bonded_to_this_carbon = [
                distanced_atom for distanced_atom in distanced_nonhydrogen_list[1:5]
                if self.measure_atom_atom_dist(int(carbon_hash), distanced_atom['#']) < self.bond_deciding_distance]
            print("non-H atoms bonded to atom %s: %s" % (
                carbon_hash, str([nonhydrogen['#'] for nonhydrogen in nonhydrogens_bonded_to_this_carbon])
            ))

            for nonhydrogen_bonded_to_this_carbon in nonhydrogens_bonded_to_this_carbon:
                # If this carbon has a neighbour not in the sp2 1e list, we know this must be a 2e carbon!
                if str(nonhydrogen_bonded_to_this_carbon['#']) not in sp2_1e_replacement_list:
                    if str(carbon_hash) not in sp2_2e_replacement_list:
                        sp2_2e_replacement_list.append(str(carbon_hash))

        # Now go back and remove 2e pseudocarbons from the 1e list.
        for carbon_hash in sp2_2e_replacement_list:
                if str(carbon_hash) in sp2_1e_replacement_list:
                    sp2_1e_replacement_list.remove(str(carbon_hash))

        log_file = open(os.path.join(os.path.dirname(self.coord_file_path), 'pseudification.log'), 'w+')
        log_file.writelines(
            '1e sp2 carbon indices: %s \n2e sp2 carbon indices: %s \nsp3 carbon indices: %s \nsp3 2e carbon indices: %s\n' % (
                ','.join(sp2_hash for sp2_hash in sp2_1e_replacement_list),
                ','.join(sp2_hash for sp2_hash in sp2_2e_replacement_list),
                ','.join(sp3_hash for sp3_hash in sp3_replacement_list),
                ','.join(sp3_hash for sp3_hash in sp3_2e_replacement_list)
                ))

        sp2_1e_coord_command = 'mn sp2 %s' % (','.join(sp2_1e_replacement_list))
        print("sp2 command: %s" % sp2_1e_coord_command)
        sp2_2e_coord_command = 'mn sp2 %s 2e' % (','.join(sp2_2e_replacement_list))
        print("sp2 command: %s" % sp2_2e_coord_command)
        sp3_coord_command = 'mn sp3 %s' % (','.join(sp3_replacement_list))
        print("sp3 command: %s" % sp3_coord_command)
        sp3_2e_coord_command = 'mn sp3 %s' % (','.join(sp3_2e_replacement_list))
        print("sp3 2e command: %s" % sp3_2e_coord_command)

        if 'nosp3' not in sysargs:
            self.pseudopotentialise_ethane_like_molecule(sp3_coord_command.split(), execute_deletion=False)
            self.pseudopotentialise_ethane_like_molecule(sp3_2e_coord_command.split(), execute_deletion=False)
        if sp2_1e_replacement_list:
            self.pseudopotentialise_molecule(sp2_1e_coord_command.split(), execute_deletion=False)
        if sp2_2e_replacement_list:
            self.pseudopotentialise_molecule(sp2_2e_coord_command.split(), execute_deletion=False)

        self.delete_specified_atoms(sp2_deletion_list + sp3_deletion_list + sp3_2e_deletion_list)

        print("Identifying 2-electron sp2 carbons...")
        # Now need to work out where the 2e sp2 carbons are
        self.coord_list = []
        self.read_coords()
        carbon_atoms = [atom for atom in self.coord_list if atom["el"] == 'c']
        sp2_pseudocarbon_list = []

        for atom in carbon_atoms:
            carbon_pseudos = self.identify_pseudocarbon_potentials(atom['#'])
            # if 6 atoms within pseudo-distance this is an sp2 pseudo-carbon
            if len(carbon_pseudos) == 6:
                sp2_pseudocarbon_list.append(atom)
        print("Re-discovered %s sp2 carbons." % str(len(sp2_pseudocarbon_list)))

        # Now check for ncore=4 sp2 pseudocarbons
        pseudopotential_hashes_to_delete = []
        for atom in sp2_pseudocarbon_list:
            distanced_nonpsnonh_list = self.order_atoms_by_distance_from(atom['#'], element='!pseudo')
            nonpsnonhs_bonded_to_this_atom = [distanced_atom for distanced_atom in distanced_nonpsnonh_list[1:5] if
                                                self.measure_atom_atom_dist(atom['#'], distanced_atom['#']) < self.bond_deciding_distance]
            print("Carbons bonded to atom %s: %s" % (str(atom['#']),
                                                     str([nonpsnonh['#'] for nonpsnonh in nonpsnonhs_bonded_to_this_atom])))

            for nonpsnonh_bonded_to_this_atom in nonpsnonhs_bonded_to_this_atom:
                # If this carbon has a neighbour not in the sp2 list, we know this must be a 2e carbon!
                if nonpsnonh_bonded_to_this_atom not in sp2_pseudocarbon_list:
                    def distance_from(list_atom):
                        return self.measure_atom_atom_dist(nonpsnonh_bonded_to_this_atom['#'], list_atom['#'])
                    carbon_pseudos = self.identify_pseudocarbon_potentials(atom['#'])
                    # find pseudos closest to the other atom
                    pseudos_distanced_from_sp2_2e = sorted(carbon_pseudos, key=distance_from)
                    pseudopotential_hashes_to_delete.append(pseudos_distanced_from_sp2_2e[0]['#'])
                    pseudopotential_hashes_to_delete.append(pseudos_distanced_from_sp2_2e[1]['#'])

        self.delete_specified_atoms(pseudopotential_hashes_to_delete)

        # Read final coordinates
        self.coord_list = []
        self.read_coords()
        carbon_atoms = [atom for atom in self.coord_list if atom["el"] == 'c']
        sp2_pseudocarbon_list = []
        sp2_2e_pseudocarbon_list = []
        sp2_2e_pseudohydrogen_list = []
        sp3_pseudocarbon_list = []
        sp3_2e_pseudocarbon_list = []

        for atom in carbon_atoms:
            carbon_pseudos = self.identify_pseudocarbon_potentials(atom['#'])

            # if 2 atoms within pseudo-distance this is an sp3 2e pseudo-carbon
            if len(carbon_pseudos) == 2:
                sp3_2e_pseudocarbon_list.append(atom)

            # if 3 atoms within pseudo-distance this is an sp3 pseudo-carbon
            if len(carbon_pseudos) == 3:
                sp3_pseudocarbon_list.append(atom)

            # if 4 atoms within pseudo-distance this is an sp2 2e pseudo-carbon
            elif len(carbon_pseudos) == 4:
                sp2_2e_pseudocarbon_list.append(atom)
                sp2_2e_pseudohydrogen_list.extend(carbon_pseudos)

            # if 6 atoms within pseudo-distance this is an sp2 pseudo-carbon
            elif len(carbon_pseudos) == 6:
                sp2_pseudocarbon_list.append(atom)


        log_file.writelines(
            'sp2 pseudocarbon indices: %s \nsp3 pseudocarbon indices: %s\nsp3 2e pseudocarbon indices: %s\nsp2 2e pseudocarbon indices: %s\nsp2 2e pseudohydrogen indices: %s\n' % (
                ','.join(str(carbon['#']) for carbon in sp2_pseudocarbon_list),
                ','.join(str(carbon['#']) for carbon in sp3_pseudocarbon_list),
                ','.join(str(carbon['#']) for carbon in sp3_2e_pseudocarbon_list),
                ','.join(str(carbon['#']) for carbon in sp2_2e_pseudocarbon_list),
                ','.join(str(carbon['#']) for carbon in sp2_2e_pseudohydrogen_list)
                ))

        # Need to supply potentials to atoms
        define_cmds_path = os.path.join(os.path.dirname(self.coord_file_path), 'define_add_pseudos')
        with open(os.path.join(define_cmds_path), 'w') as var_file:
            var_file.writelines(define_cmds % (
                # all atoms
                self.supply_ecps_bases_to_define('all', 'b', self.all_e_basis),
                # sp2 potentials
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp2_pseudocarbon_list], 'b', self.sp2_pseudo_carbon_basis),
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp2_pseudocarbon_list], 'ecp', self.sp2_carbon_ecp),
                self.supply_ecps_bases_to_define(self.sp2_pseudo_element, 'b', 'none'),
                self.supply_ecps_bases_to_define(self.sp2_pseudo_element, 'ecp', self.sp2_hydrogen_ecp),
                # sp3 potentials
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp3_pseudocarbon_list], 'b', self.sp3_pseudo_carbon_basis),
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp3_pseudocarbon_list], 'ecp', self.sp3_carbon_ecp),
                self.supply_ecps_bases_to_define(self.sp3_pseudo_element, 'b', 'none'),
                self.supply_ecps_bases_to_define(self.sp3_pseudo_element, 'ecp', self.sp3_hydrogen_ecp),
                # sp2 2e potentials
                self.supply_ecps_bases_to_define(self.sp2_pseudo_element, 'b', 'none'),
                self.supply_ecps_bases_to_define([hydrogen['#'] for hydrogen in sp2_2e_pseudohydrogen_list], 'ecp', self.sp2_2e_hydrogen_ecp),
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp2_2e_pseudocarbon_list], 'b', self.sp2_pseudo_carbon_basis),
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp2_2e_pseudocarbon_list], 'ecp', self.sp2_2e_carbon_ecp),
                # sp3 2e potentials
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp3_2e_pseudocarbon_list], 'b', self.sp3_2e_pseudo_carbon_basis),
                self.supply_ecps_bases_to_define([carbon['#'] for carbon in sp3_2e_pseudocarbon_list], 'ecp', self.sp3_2e_carbon_ecp),
                self.supply_ecps_bases_to_define(self.sp3_2e_pseudo_element, 'b', 'none'),
                self.supply_ecps_bases_to_define(self.sp3_2e_pseudo_element, 'ecp', self.sp3_2e_hydrogen_ecp)
            ))

        self.run_define('define_add_pseudos')

    def pseudopotentialise_molecule(self, sysargs=None, execute_deletion=True):
        """Creates sets of pseudo-potentials for specified C atoms a la CH3 radical."""

        # Find atoms to replace
        deletion_list = []
        if len(sysargs) > 2:
            if 'del' in sysargs:
                deletion_list = self.parse_coord_list(sysargs[4])
            replacement_list = self.parse_coord_list(sysargs[2])
            atoms_to_potentialise = list(item for item in self.coord_list if item["#"] in replacement_list)
        else:
            atoms_to_potentialise = (item for item in self.coord_list if item["el"] == 'c')
            deletion_list = (item for item in self.coord_list if item["el"] == 'h')
        print('Pseudo-potentialising carbon atoms %s ...' % [atom['#'] for atom in atoms_to_potentialise])

        potential_coords_list = []

        for atom in atoms_to_potentialise:
            distanced_atom_list = self.order_atoms_by_distance_from(atom['#'])
            distanced_carbon_list = self.order_atoms_by_distance_from(atom['#'], element='c')

            atoms_within_bond_distance = self.get_atoms_within_bond_distance(atom['#'])

            if len(distanced_carbon_list) in (1, 2):
                primary_vector = None
                for non_c_atom in atoms_within_bond_distance:
                    if non_c_atom['el'] != 'h':
                        primary_vector = self.vectorise_atom(non_c_atom['#']) - self.vectorise_atom(atom['#'])
                if primary_vector is None:
                    primary_vector = self.vectorise_atom(distanced_atom_list[1]['#']) - self.vectorise_atom(atom['#'])

                normal_vector = numpy.cross(
                    self.vectorise_atom(distanced_atom_list[1]['#']) - self.vectorise_atom(atom['#']),
                    self.vectorise_atom(distanced_atom_list[2]['#']) - self.vectorise_atom(atom['#'])
                )
            else:
                primary_vector = self.vectorise_atom(distanced_carbon_list[1]['#']) - self.vectorise_atom(atom['#'])

                normal_vector = numpy.cross(
                    self.vectorise_atom(distanced_carbon_list[1]['#']) - self.vectorise_atom(atom['#']),
                    self.vectorise_atom(distanced_carbon_list[2]['#']) - self.vectorise_atom(atom['#'])
                )

            if '2e' in sysargs:
                print('2e keyword used. Using 2e-sp2 parameters.')
                primary_potential_vector = self.lengtherise_vector(primary_vector,
                                                                   self.sp2_2e_atom_potential_set_distance)
                potential_set_split_vector = self.lengtherise_vector(normal_vector,
                                                                     self.sp2_2e_potential_set_split_distance)
            else:
                primary_potential_vector = self.lengtherise_vector(primary_vector,
                                                                   self.sp2_1e_atom_potential_set_distance)
                potential_set_split_vector = self.lengtherise_vector(normal_vector,
                                                                     self.sp2_1e_potential_set_split_distance)

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

            if self.add_primary_vector_potentials_as_coords is False:
                del relative_potential_vectors[0]
                del relative_potential_vectors[0]

            # potential coords are still relative to their atom, now make them real.
            for vector in relative_potential_vectors:
                potential_coords_list.append(
                    {'#': 0, 'el': self.sp2_pseudo_element, 'x': vector[0]+atom['x'], 'y': vector[1]+atom['y'], 'z': vector[2]+atom['z']},
                )

        # Now add potentials to coord list, after removing the 'real' hydrogen atoms.
        if execute_deletion is True:
            self.delete_specified_atoms(deletion_list)
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=False)

    def pseudopotentialise_ethane_like_molecule(self, sysargs, execute_deletion=True):
        """Creates pseudo-potentials to replace specified atoms a la ethane."""

        # Find atoms to replace
        deletion_list = []
        potential_coords_list = []
        if len(sysargs) > 2:
            if 'del' in sysargs:
                deletion_list = self.parse_coord_list(sysargs[4])
            replacement_list = self.parse_coord_list(sysargs[2])
            atoms_to_replace = list(item for item in self.coord_list if item["#"] in replacement_list)
        else:
            atoms_to_replace = (item for item in self.coord_list if item["el"] == 'c')
            deletion_list = (item for item in self.coord_list if item["el"] == 'h')
        print('Pseudo-potentialising atoms %s ...' % [atom['#'] for atom in atoms_to_replace])

        # Option to place a potential on the *opposite* side of the carbon as well.
        dipolar_potentials = False
        if 'dipole' in sysargs:
            print('Dipolar potentialisation activated...')
            dipolar_potentials = True

        for atom in atoms_to_replace:
            # Find vector from nearest carbon.
            distanced_carbon_list = self.order_atoms_by_distance_from(atom['#'], element='c')

            vector_from_nearest_carbon = self.vectorise_atom(atom['#']) \
                - self.vectorise_atom(distanced_carbon_list[0]['#'])
            vector_to_nearest_carbon = self.vectorise_atom(distanced_carbon_list[0]['#']) \
                - self.vectorise_atom(atom['#'])

            # Lengtherise vector from carbon to give relative pp coordinates.
            vector_c_to_new_pp = self.lengtherise_vector(vector_from_nearest_carbon, self.sp3_1e_atom_potential_distance)
            vector_c_to_new_dipole_pp = self.lengtherise_vector(vector_to_nearest_carbon, self.sp3_1e_atom_potential_distance)

            # Add to carbon coords to get new pp coords.
            potential_coords_list.append(
                {'#': 0, 'el': self.sp3_pseudo_element,
                 'x': vector_c_to_new_pp[0] + distanced_carbon_list[0]['x'],
                 'y': vector_c_to_new_pp[1] + distanced_carbon_list[0]['y'],
                 'z': vector_c_to_new_pp[2] + distanced_carbon_list[0]['z']},
            )
            if dipolar_potentials is True:
                # Add to carbon coords to get new pp coords.
                potential_coords_list.append(
                    {'#': 0, 'el': self.sp3_pseudo_element,
                     'x': vector_c_to_new_dipole_pp[0] + distanced_carbon_list[0]['x'],
                     'y': vector_c_to_new_dipole_pp[1] + distanced_carbon_list[0]['y'],
                     'z': vector_c_to_new_dipole_pp[2] + distanced_carbon_list[0]['z']},
                )

        # Now add potentials to coord list, after removing the 'real' atoms.
        if execute_deletion is True:
            self.delete_specified_atoms(deletion_list)
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=False)

    def pseudopotentialise_sp3_carbon(self, sysargs, execute_deletion=True):
        """Potentialises an sp3 atom by reference to its carbon, unlike the other ethane function."""
        #TODO
        # Find atoms to replace
        deletion_list = []
        potential_coords_list = []
        if len(sysargs) > 2:
            if 'del' in sysargs:
                deletion_list = self.parse_coord_list(sysargs[4])
            replacement_list = self.parse_coord_list(sysargs[2])
            atoms_to_replace = list(item for item in self.coord_list if item["#"] in replacement_list)
            atoms_to_potentialise = self.parse_coord_list(sysargs[2])
        else:
            atoms_to_potentialise = (item for item in self.coord_list if item["el"] == 'c')
            #deletion_list = (item for item in self.coord_list if item["el"] == 'h')
        print('Pseudo-potentialising atoms %s ...' % [atom['#'] for atom in atoms_to_potentialise])

        for atom in atoms_to_potentialise:
            three_nearest_hydrogens = self.order_atoms_by_distance_from(atom['#'], element='h')


            # Find vector from nearest carbon.
            distanced_carbon_list = self.order_atoms_by_distance_from(atom['#'], element='c')

            vector_from_nearest_carbon = self.vectorise_atom(atom['#']) \
                - self.vectorise_atom(distanced_carbon_list[0]['#'])
            vector_to_nearest_carbon = self.vectorise_atom(distanced_carbon_list[0]['#']) \
                - self.vectorise_atom(atom['#'])

            # Lengtherise vector from carbon to give relative pp coordinates.
            vector_c_to_new_pp = self.lengtherise_vector(vector_from_nearest_carbon, self.sp3_1e_atom_potential_distance)

            # Add to carbon coords to get new pp coords.
            potential_coords_list.append(
                {'#': 0, 'el': self.sp3_pseudo_element,
                 'x': vector_c_to_new_pp[0] + distanced_carbon_list[0]['x'],
                 'y': vector_c_to_new_pp[1] + distanced_carbon_list[0]['y'],
                 'z': vector_c_to_new_pp[2] + distanced_carbon_list[0]['z']},
            )

        # Now add potentials to coord list, after removing the 'real' atoms.
        if execute_deletion is True:
            self.delete_specified_atoms(deletion_list)
        for potential_coord in potential_coords_list:
            self.write_coord(potential_coord, overwrite=False)


if __name__ == "__main__":
    control = CoordControl()
    control.read_coords()
    if sys.argv[1] == 'sp2':
        control.pseudopotentialise_molecule(sys.argv)
    elif sys.argv[1] == 'sp3':
        control.pseudopotentialise_ethane_like_molecule(sys.argv)
    elif sys.argv[1] == 'guess':
        control.guess_potentialisation(sys.argv)
    elif sys.argv[1] == 'repotentialise_sp2':
        for pseudo_atom in control.parse_coord_list(sys.argv[2]):
            control.repseudopotentialise_sp2_atom(pseudo_atom, float(sys.argv[3]), float(sys.argv[4]))

        print('Re-potentialised atom %s, with set distance %s, split %s' % (control.parse_coord_list(sys.argv[2]),
                                                                            float(sys.argv[3]),
                                                                            float(sys.argv[4])))
    elif sys.argv[1] == 'repotentialise_sp3':
        for pseudo_atom in control.parse_coord_list(sys.argv[2]):
            control.set_potential_distance_to(pseudo_atom, float(sys.argv[3]))

        print('Re-potentialised atoms %s, with set distance %s' % (control.parse_coord_list(sys.argv[2]),
                                                                   float(sys.argv[3])))

    else:
        print('Incorrect sysargs given.')


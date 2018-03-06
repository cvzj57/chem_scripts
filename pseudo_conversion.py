import subprocess
import sys
import shutil
from chem_lib.optimise import BasisControl
from chem_lib.coord import CoordControl
import os

define_cmds = '''

a coord
ired
*
b "he" none
ecp
%s
ecp-p
file basis
ecp "he" ecp-s
file basis
*
q

dft
func pbe0
grid m4
on
q
q
'''

to_tddft_cmds = '''






ex
cist
*
a 1
q
q


*
'''


class ConversionControl:
    def __init__(self):
        self.all_electron_calculation_folder_path = ''
        self.pseudo_calculation_folder_path = ''
        self.pseudo_basis_path = ''
        self.pseudo_carbon_indices = []
        self.pseudo_hydrogen_indices = []
        self.initial_homo_index = 0
        self.current_homo_indices = []
        self.basis = BasisControl()
        self.hybridisation = 'sp2'
        self.folder_names = [
            'singlet',
            'triplet',
            '1sti',
            'tddft',
        ]

    def getopts(argv):
        opts = {}  # Empty dictionary to store key-value pairs.
        while argv:  # While there are arguments left to parse...
            if argv[0][0] == '-':  # Found a "-name value" pair.
                opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
            argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
        return opts

    def run_define(self, cmds):
        command = 'define < %s' % cmds
        print(command)
        subprocess.call(command, shell=True, cwd=self.pseudo_calculation_folder_path)

    def delete_occupations_control(self):
        occupation_keys = ['scfmo', 'uhfmo_alpha', 'uhfmo_beta', 'uhf', 'alpha shells', 'beta shells', 'closed shells',
                           'scfiterlimit']
        for occupation_key in occupation_keys:
            subprocess.call('kdg %s' % occupation_key, cwd=self.pseudo_calculation_folder_path, shell=True)

    def modify_atoms_control(self):

        define_cmds_path = os.path.join(self.pseudo_calculation_folder_path, 'define_add_pseudos')
        shutil.copy(os.path.join(self.pseudo_basis_path, 'basis'), self.pseudo_calculation_folder_path)
        os.remove(os.path.join(self.pseudo_calculation_folder_path, 'control'))

        with open(os.path.join(define_cmds_path), 'w') as var_file:
            var_file.writelines(define_cmds % ','.join(self.pseudo_carbon_indices))

        self.run_define('define_add_pseudos')

    def add_excitation_control(self):

        to_tddft_cmds_path = os.path.join(self.pseudo_calculation_folder_path, 'define_add_excitation')

        with open(os.path.join(to_tddft_cmds_path), 'w') as var_file:
            var_file.writelines(to_tddft_cmds)

        self.run_define('define_add_excitation')

    def add_occupation_to_control(self, homo_indices):

        var_file = open(os.path.join(self.pseudo_calculation_folder_path, 'control'), 'r')
        var_file_data = var_file.readlines()

        if len(homo_indices) > 1:
            var_file_data.insert(
                -1,
                "$uhf\n$uhfmo_alpha none file=alpha\n$uhfmo_beta none file=beta\n"
                "$alpha shells\n a       1-%s         ( 1 )\n"
                "$beta shells\n a       1-%s         ( 1 )\n"
                "$scfiterlimit 300\n"
                % (homo_indices[0], homo_indices[1])
            )
        else:
            var_file_data.insert(
                -1,
                "$scfmo none file=mos \n$closed shells \n  a       1-%s    ( 2 )\n"
                "$scfiterlimit 300\n" % homo_indices[0]
            )

        # And write everything back
        with open(os.path.join(self.pseudo_calculation_folder_path, 'control'), 'w') as var_file:
            if var_file_data:
                var_file.writelines(var_file_data)
        var_file.close()

    def perform_conversion(self, folder_paths, run_escf=False):

        self.all_electron_calculation_folder_path = folder_paths[0]
        self.pseudo_calculation_folder_path = folder_paths[1]

        print("converting %s to %s..." % (self.all_electron_calculation_folder_path,
                                          self.pseudo_calculation_folder_path))

        # copy everything to new directory
        try:
            shutil.copytree(self.all_electron_calculation_folder_path, self.pseudo_calculation_folder_path)
        except Exception as e:
            print('What are you playing at?')

        print('files copied...')
        print(self.pseudo_hydrogen_indices)

        coord_control = CoordControl()
        coord_control.coord_file_path = os.path.join(self.pseudo_calculation_folder_path, 'coord')
        coord_control.read_coords()

        if self.hybridisation == 'sp2':
            coord_command = 'mn %s %s del %s' % (self.hybridisation,
                                                 ' '.join(self.pseudo_carbon_indices),
                                                 ' '.join(self.pseudo_hydrogen_indices))
            print('coord command: %s' % coord_command)
            coord_control.pseudopotentialise_molecule(coord_command.split())
        elif self.hybridisation == 'sp3':
            coord_command = 'mn %s %s del %s' % (self.hybridisation,
                                                 ' '.join(self.pseudo_hydrogen_indices),
                                                 ' '.join(self.pseudo_hydrogen_indices))
            print('coord command: %s' % coord_command)
            coord_control.pseudopotentialise_ethane_like_molecule(coord_command.split())

        # run define and add new ecps/bases
        self.modify_atoms_control()
        # remove current orbital occupations
        self.delete_occupations_control()
        # add new occupations to control
        self.add_occupation_to_control(self.current_homo_indices)

        # run new calculation
        if run_escf is True:
            self.basis.run_dscf(add_to_log=True, file_path=self.pseudo_calculation_folder_path)
            os.remove(os.path.join(self.pseudo_calculation_folder_path, 'cist_a'))
            self.add_excitation_control()
            self.basis.run_dscf(add_to_log=True, file_path=self.pseudo_calculation_folder_path)
            self.basis.run_escf(add_to_log=True, file_path=self.pseudo_calculation_folder_path)
        else:
            self.basis.run_dscf(add_to_log=True, file_path=self.pseudo_calculation_folder_path)

    def run(self, sysargs):
        reference_calcs_path = sysargs[1]
        dest_calcs_path = sysargs[2]

        self.pseudo_hydrogen_indices = input('Enter pseudo-hydrogen indices: ').split(' ')

        self.pseudo_carbon_indices = input('Enter pseudo-carbon indices: ').split(' ')

        self.initial_homo_index = int(input('Enter HOMO index: '))

        # run singlet
        self.current_homo_indices = [self.initial_homo_index]
        self.perform_conversion([os.path.join(reference_calcs_path, 'singlet'),
                                 os.path.join(dest_calcs_path, 'singlet')])
        # run triplet
        self.current_homo_indices = [self.initial_homo_index + 1, self.initial_homo_index - 1]
        self.perform_conversion([os.path.join(reference_calcs_path, 'triplet'),
                                 os.path.join(dest_calcs_path, 'triplet')])
        # run 1sti
        self.current_homo_indices = [self.initial_homo_index, self.initial_homo_index - 1]
        self.perform_conversion([os.path.join(reference_calcs_path, '1sti'),
                                 os.path.join(dest_calcs_path, '1sti')])
        # run tddft
        self.current_homo_indices = [self.initial_homo_index]
        self.perform_conversion([os.path.join(reference_calcs_path, 'tddft'),
                                 os.path.join(dest_calcs_path, 'tddft')],
                                run_escf=True)


if __name__ == "__main__":
    control = ConversionControl()
    control.run(sys.argv)

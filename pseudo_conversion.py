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


class ConversionControl:
    def __init__(self):
        self.all_electron_calculation_folder_path = ''
        self.pseudo_calculation_folder_path = ''
        self.pseudo_basis_path = ''
        self.pseudo_carbon_indices = []
        self.pseudo_hydrogen_indices = []
        self.basis = BasisControl()
        self.hybridisation = 'sp3'

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

    def run(self, sysargs):

        self.all_electron_calculation_folder_path = sysargs[1]
        self.pseudo_calculation_folder_path = sysargs[2]

        # copy everything to new directory
        try:
            shutil.copytree(self.all_electron_calculation_folder_path, self.pseudo_calculation_folder_path)
        except Exception as e:
            print('What are you playing at?')

        print('files copied...')
        coord_command = input('Enter pseudo-hydrogen indices: ')
        self.pseudo_hydrogen_indices = ['mn', self.hybridisation] + coord_command.split(' ') + ['del'] + coord_command.split(' ')
        print(self.pseudo_hydrogen_indices)

        coord_control = CoordControl()
        coord_control.coord_file_path = os.path.join(self.pseudo_calculation_folder_path, 'coord')
        coord_control.read_coords()

        if self.hybridisation == 'sp2':
            coord_control.pseudopotentialise_molecule()
        elif self.hybridisation == 'sp3':
            coord_control.pseudopotentialise_ethane_like_molecule(self.pseudo_hydrogen_indices)

        # run define and add new ecps/bases
        self.pseudo_carbon_indices = input('Enter pseudo-carbon indices: ').split(' ')
        self.modify_atoms_control()
        # remove current orbital occupations
        self.delete_occupations_control()
        # add new occupations to control
        homo_indices = input('Enter HOMO indices: ').split(' ')
        self.add_occupation_to_control(homo_indices)

        # run new calculation
        self.basis.run_dscf(add_to_log=True, file_path=self.pseudo_calculation_folder_path)


if __name__ == "__main__":
    control = ConversionControl()
    control.run(sys.argv)

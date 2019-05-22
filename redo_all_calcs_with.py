'''Re-running everything with new potentials/functionals/etc is a pain. Bash scripting of only some help.
This script should, with the correct parameters, be able to re-run all the CnHn and CnHn+2 calculations with whatever
parameters I want.

I'll start with the basic folder structure with the geometries already optimised at the UHF singlet level.

It might be nice to have a few steps in the process that can be turned on and off, or maybe different functions
for different kinds of modifications e.g. a new functional or pseudo-potential.'''

import os
import shutil
from chem_lib.coord import CoordControl
import subprocess
import sys
import json

example_info_dict = {
    'coord_file_name': 'c2h4',
    'atoms_to_potentialise': '1,2',
    'remaining_electrons': 2
}

ch_sym_cmds = '''
y
sy %s
*
no




*






*
*
'''

#TODO this means i will need the carbons first in all coord files?

class CalculationRunner:
    def __init__(self):
        self.path_to_cnhn_2_calcs = '.'
        self.path_to_cnhn_calcs = '.'
        self.pp_basis_file = 'test_basis'
        self.functional_list = ['pbe', 'pbe0', 'tpss', 'tpssh', 'b3-lyp']
        self.occupation_setup_list = ['singlet_uhf', 'triplet', 'cation', 'tddft']
        self.potentialisation = 'alpha'
        self.symmetry = 'c1'
        self.run_jobs = False
        self.qsub_command = 'qsub submit.job'
        self.coordlib_path = 'coord_lib'
        with open(os.path.dirname(os.path.realpath(__file__))+'/chem_lib/calc_molecular_library.json', 'r') as json_file:
            self.molecule_library = json.load(json_file)


    def set_symmetry(self, calc_path, new_symmetry):
        # feed script into turbomole with new symmetry
        command = 'define < %s' % (ch_sym_cmds % new_symmetry)
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=calc_path)
        return

    def add_occupation(self, calc_path, total_es, occ_type, orbital_symgroup='a'):
        control_file = open(os.path.join(calc_path, 'control'), 'w+')
        control_file_data = control_file.readlines()

        alpha = 1
        beta = 1

        if occ_type == 'singlet_uhf':
            alpha, beta=total_es/2
        elif occ_type == 'triplet':
            alpha=total_es/2+1
            beta=total_es/2-1
        elif occ_type == 'cation':
            alpha=total_es/2
            beta=total_es/2-1

        control_file_data.insert(-1, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % ('$uhf',
                                                                       '$uhfmo_alpha none file=alpha',
                                                                       '$uhfmo_beta none file=beta',
                                                                       '$alpha shells',
                                                                       ' %s %s (1)' % (orbital_symgroup,
                                                                                       '1-%s' % str(alpha)),
                                                                       '$beta shells',
                                                                       ' %s %s (1)' % (orbital_symgroup,
                                                                                       '1-%s' % str(beta))))
        return

    def add_dft_functional(self, calc_path, functional_name):
        control_file = open(os.path.join(calc_path, 'control'), 'w+')
        control_file_data = control_file.readlines()
        control_file_data.insert(-1, '%s\n%s\n%s\n' % ('$dft', ' func %s' % functional_name, ' grid m4'))
        control_file.writelines(control_file_data)
        control_file.close()
        return

    def run(self, supplied_basis_file, molecule_set, set_dist, split_dist=0.25):
        # read list of molecules to create
        molecule_list = self.molecule_library[molecule_set]
        molecule_list = example_info_dict

        # read potentials to use

        # create folder structure with s,t etc and all functionals and copy all-e coords into folders

        coord = CoordControl()

        for molecule in molecule_list:
            molname = molecule['coord_file_name']
            this_script_path = os.path.realpath(__file__)
            os.mkdir(molname)
            for functional_name in self.functional_list:
                singlet_path = os.path.join(molname, functional_name, 'singlet_uhf')

                os.mkdir(os.path.join(molname, functional_name))
                os.mkdir(os.path.join(molname, functional_name, 'singlet_uhf'))
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, molname), singlet_path)
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, 'submit.job'), singlet_path)
                shutil.copy(supplied_basis_file, singlet_path)

                # run the guess script for each coord file
                parsed_atoms_to_potentialise = coord.parse_coord_list(molecule['atoms_to_potentialise'])
                coord.guess_potentialisation(['elephants', 'guess', 'incl', parsed_atoms_to_potentialise])

                # repotentialise with correct values
                if self.potentialisation == 'alpha':
                    for pseudo_atom in parsed_atoms_to_potentialise:
                        coord.repseudopotentialise_sp2_atom(pseudo_atom, set_dist, split_dist)
                elif self.potentialisation == 'gamma':
                    for pseudo_atom in parsed_atoms_to_potentialise:
                        coord.set_potential_distance_to(pseudo_atom, set_dist)
                # set pseudopotential values
                shutil.copy(supplied_basis_file, singlet_path)

                # set any extra parameters functional
                self.add_dft_functional(os.path.join(molname, functional_name, 'singlet_uhf'), functional_name)

                # set symmetry
                self.set_symmetry(os.path.join(molname, functional_name, 'singlet_uhf'), self.symmetry)

                # copy to triplet, cation, tddft and set occupations
                triplet_path = os.path.join(molname, functional_name, 'triplet')
                cation_path = os.path.join(molname, functional_name, 'cation')
                tddft_path = os.path.join(molname, functional_name, 'tddft')

                if 'triplet' in self.occupation_setup_list:
                    shutil.copytree(singlet_path, triplet_path)
                    self.add_occupation(triplet_path, molecule['total_es'], 'triplet', 'a')

                if 'cation' in self.occupation_setup_list:
                    shutil.copytree(singlet_path, cation_path)
                    self.add_occupation(cation_path, molecule['total_es'], 'cation', 'a')

                if 'tddft' in self.occupation_setup_list:
                    shutil.copytree(singlet_path, tddft_path)
                    self.add_occupation(tddft_path, molecule['total_es'], 'singlet_uhf', 'a')

                # Now the singlet occ
                self.add_occupation(singlet_path, molecule['total_es'], 'singlet_uhf', 'a')

                # sort out tddft escf stuff
                if 'tddft' in self.occupation_setup_list:
                    pass
                #TODO don't forget submit.job needs escf

                # run calcs if desired
                if self.run_jobs is True:
                    subprocess.call(self.qsub_command, shell=True, cwd=singlet_path)
                    if 'triplet' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=triplet_path)
                    if 'triplet' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=cation_path)
                    if 'triplet' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=tddft_path)

        # gather results. what's the best format?
        # grep from dscf?
        #os.path.basename(your_path)
        # parse and insert into csv? Can you copy macros straight from python into csv?
        # needs to end up in latex table format
        # If I can read it sensibly enough into python I could do the averaging and percentage errors manually

        # subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return

if __name__ == "__main__":
    calcrunner = CalculationRunner()
    args = sys.argv.split()
    calcrunner.run(args)
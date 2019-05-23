'''Re-running everything with new potentials/functionals/etc is a pain. Bash scripting of only some help.
This script should, with the correct parameters, be able to re-run all test set calculations with whatever
parameters I want.

Uses libraries of coords and specifications of potential occupations to build geometries, add whatever
potentials and run.

Even gathers results into a handy CSV!
'''

import os
import shutil
from chem_lib.coord import CoordControl
import subprocess
import sys
import json

example_mol_list = [{
    'coord_file_name': 'c2h4',
    'atoms_to_potentialise': '1,2',
    'total_es': 2
}]

ch_sym_cmds = '''
y
sy %s
*
no




*






*
*
'''

ref_to_t_cmds = '''


y
eht


n
t
*




*
'''

ref_to_i_cmds = '''



y
eht

1






*
'''

setup_cmds = '''


a coord
*
no
*
eht


n
s
*

*
'''

#TODO this means i will need the carbons first in all coord files?

class CalculationRunner:
    def __init__(self):
        self.functional_list = ['hf', 'pbe', 'pbe0', 'tpss', 'tpssh', 'b3-lyp']
        self.occupation_setup_list = ['singlet_uhf', 'triplet', 'cation', 'tddft']
        self.potentialisation = 'alpha'
        self.symmetry = 'cs'
        self.run_jobs = True
        self.qsub_command = 'qsub submit.job'
        self.coordlib_path = 'coord_lib'
        self.occ_orb_sym = 'a\"'
        self.ex_symmetry = 'a\''
        self.scfinstab = 'urpa'
        with open(os.path.dirname(os.path.realpath(__file__))+'/coord_lib/molecular_library.json', 'r') as json_file:
            self.molecule_library = json.load(json_file)

    def sort_atoms(self, coord_file_name):
        coord_file = open(coord_file_name, 'r')
        coord_file_data = coord_file.readlines()
        h_lines = []
        for line in coord_file_data:
            if ' h' in line:
                h_lines.append(line)

        for line in h_lines:
            coord_file_data.insert(-2, line)

        for line in h_lines:
            coord_file_data.remove(line)

        with open(coord_file_name, 'w') as coord_file:
            coord_file.writelines(coord_file_data)

        print('Moved hydrogens to end of coord file.')
        coord_file.close()
        return

    @staticmethod
    def increase_scfiterlimit(calc_path):
        command = 'kdg scfiterlimit'
        subprocess.call(command, shell=True, cwd=calc_path)
        command = 'adg scfiterlimit 3000'
        subprocess.call(command, shell=True, cwd=calc_path)

    def add_tddft_excitation(self, calc_path, ex_symmetry, instab='rpas'):
        control_file = open(os.path.join(calc_path, 'control'), 'r')
        control_file_data = control_file.readlines()
        control_file_data.insert(-1, '%s\n%s\n%s\n%s\n%s\n' % ('$scfinstab  %s' % instab,
                                                               '$soes',
                                                               ' %s 1' % ex_symmetry,
                                                               '$rpacor 1000',
                                                               '$escfiterlimit 3000'))
        with open(os.path.join(calc_path, 'control'), 'w') as control_file:
            control_file.writelines(control_file_data)
        control_file.close()

        submit_file = open(os.path.join(calc_path, 'submit.job'), 'r')
        submit_file_data = submit_file.readlines()
        submit_file_data.insert(-5, 'escf > escf.out\n')
        with open(os.path.join(calc_path, 'submit.job'), 'w') as submit_file:
            submit_file.writelines(submit_file_data)
        submit_file.close()

        return

    def set_symmetry(self, calc_path, new_symmetry):
        # feed script into turbomole with new symmetry
        print("Symmetry changing to %s..." % new_symmetry)
        define_cmds_path = 'define_sym_change_%s' % new_symmetry

        with open(os.path.join(calc_path, define_cmds_path), 'w') as var_file:
            var_file.writelines(ch_sym_cmds % new_symmetry)
        var_file.close()

        command = 'define < %s' % define_cmds_path
        subprocess.call(command, shell=True, cwd=calc_path)
        return

    def add_occupation(self, calc_path, total_es, occ_type, orbital_symgroup='a'):
        control_file = open(os.path.join(calc_path, 'control'), 'r')
        control_file_data = control_file.readlines()

        alpha = 1
        beta = 1

        if occ_type == 'singlet_uhf':
            alpha = beta = total_es/2
        elif occ_type == 'triplet':
            alpha = total_es/2+1
            beta = total_es/2-1
        elif occ_type == 'cation':
            alpha = total_es/2
            beta = total_es/2-1

        alpha_orb_syntax = ' %s %s (1)' % (orbital_symgroup, '1-%s' % str(int(alpha)))
        beta_orb_syntax = ' %s %s (1)' % (orbital_symgroup, '1-%s' % str(int(beta)))
        if total_es == 2 and (occ_type == 'triplet' or occ_type == 'cation'):
            beta_orb_syntax = ''

        control_file_data.insert(-1, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n' % ('$uhf',
                                                                       '$uhfmo_alpha none file=alpha',
                                                                       '$uhfmo_beta none file=beta',
                                                                       '$alpha shells',
                                                                       alpha_orb_syntax,
                                                                       '$beta shells',
                                                                       beta_orb_syntax))
        with open(os.path.join(calc_path, 'control'), 'w') as control_file:
            control_file.writelines(control_file_data)
        control_file.close()
        return

    def add_dft_functional(self, calc_path, functional_name):
        if functional_name is not 'hf':
            control_file = open(os.path.join(calc_path, 'control'), 'r')
            control_file_data = control_file.readlines()
            control_file_data.insert(-1, '%s\n%s\n%s\n' % ('$dft', ' func %s' % functional_name, ' grid m4'))
            with open(os.path.join(calc_path, 'control'), 'w') as control_file:
                control_file.writelines(control_file_data)
            control_file.close()
        return

    def run(self, supplied_basis_file, molecule_set, set_dist, split_dist=0.25):
        # read list of molecules to create
        molecule_list = self.molecule_library[molecule_set]
        #molecule_list = example_mol_list

        # read potentials to use

        # create folder structure with s,t etc and all functionals and copy all-e coords into folders

        coord = CoordControl()

        this_script_path = os.path.dirname(os.path.realpath(__file__))
        test_folder_name = 'test_%s_%s' % (molecule_set, supplied_basis_file)
        os.mkdir(test_folder_name)
        shutil.copy(supplied_basis_file, test_folder_name)
        os.chdir(test_folder_name)

        for molecule in molecule_list:
            molname = molecule['coord_file_name']
            print('Creating pseudo-%ss..' % molname)
            os.mkdir(molname)

            for functional_name in self.functional_list:
                singlet_path = os.path.join(molname, functional_name, 'singlet_uhf')

                os.mkdir(os.path.join(molname, functional_name))
                os.mkdir(os.path.join(molname, functional_name, 'singlet_uhf'))
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, molname),
                            os.path.join(singlet_path, 'coord'))
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, 'submit.job'),
                            singlet_path)
                shutil.copy(supplied_basis_file, singlet_path)

                # run the guess script for each coord file
                print('Working out where your potentials are...')
                parsed_atoms_to_potentialise = coord.parse_coord_list(molecule['atoms_to_potentialise'])
                coord.coord_file_path = os.path.join(singlet_path, 'coord')
                coord.coord_list = []
                coord.read_coords()
                coord.guess_potentialisation(['elephants', 'guess', 'incl', molecule['atoms_to_potentialise']])

                # repotentialise with correct values
                print('Repotentialising... ')
                coord.coord_list = []
                coord.read_coords()
                if molecule_set == 'alpha':
                    for pseudo_atom in parsed_atoms_to_potentialise:
                        coord.repseudopotentialise_sp2_atom(pseudo_atom, set_dist, split_dist)
                elif molecule_set == 'gamma':
                    for pseudo_atom in parsed_atoms_to_potentialise:
                        coord.set_potential_distance_to(pseudo_atom, set_dist)
                # set pseudopotential values
                shutil.copy(supplied_basis_file, singlet_path)

                # set any extra parameters functional
                print('Adding DFT functionals...')
                self.add_dft_functional(singlet_path, functional_name)

                # set symmetry
                print('Updating symmetry...')
                self.set_symmetry(singlet_path, self.symmetry)

                # increase iterlimit
                print('Increasing iterlimit...')
                self.increase_scfiterlimit(singlet_path)

                # copy to triplet, cation, tddft and set occupations
                triplet_path = os.path.join(molname, functional_name, 'triplet')
                cation_path = os.path.join(molname, functional_name, 'cation')
                tddft_path = os.path.join(molname, functional_name, 'tddft')

                if 'triplet' in self.occupation_setup_list:
                    print('Creating triplet...')
                    shutil.copytree(singlet_path, triplet_path)
                    self.add_occupation(triplet_path, molecule['total_es'], 'triplet', self.occ_orb_sym)

                if 'cation' in self.occupation_setup_list:
                    print('Creating cation...')
                    shutil.copytree(singlet_path, cation_path)
                    self.add_occupation(cation_path, molecule['total_es'], 'cation', self.occ_orb_sym)

                if 'tddft' in self.occupation_setup_list:
                    print('Creating TDDFT calc...')
                    shutil.copytree(singlet_path, tddft_path)
                    self.add_occupation(tddft_path, molecule['total_es'], 'singlet_uhf', self.occ_orb_sym)
                    self.add_tddft_excitation(tddft_path, self.ex_symmetry, self.scfinstab)

                # Now the singlet occ
                self.add_occupation(singlet_path, molecule['total_es'], 'singlet_uhf', self.occ_orb_sym)

                # run calcs if desired
                if self.run_jobs is True:
                    subprocess.call(self.qsub_command, shell=True, cwd=singlet_path)
                    if 'triplet' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=triplet_path)
                    if 'cation' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=cation_path)
                    if 'tddft' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=tddft_path)

        os.chdir('..')
        return

    def reference_run(self, molecule_set):
        # read list of molecules to create
        molecule_list = self.molecule_library[molecule_set]
        # molecule_list = example_mol_list

        this_script_path = os.path.dirname(os.path.realpath(__file__))
        test_folder_name = 'test_%s_reference' % (molecule_set)
        os.mkdir(test_folder_name)
        os.chdir(test_folder_name)

        for molecule in molecule_list:
            molname = molecule['coord_file_name']
            print('Creating all-electron-%ss..' % molname)
            os.mkdir(molname)

            for functional_name in self.functional_list:
                singlet_path = os.path.join(molname, functional_name, 'singlet_uhf')

                os.mkdir(os.path.join(molname, functional_name))
                os.mkdir(os.path.join(molname, functional_name, 'singlet_uhf'))
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, molname),
                            os.path.join(singlet_path, 'coord'))
                shutil.copy(os.path.join(this_script_path, self.coordlib_path, molecule_set, 'submit.job'),
                            singlet_path)

                # instead of all the guessing and basis/ecp assigning, we just need a define script to take us through a standard setup
                # then a couple more for triplet/cation setups
                self.alle_singlet_uhf_setup(singlet_path)

                # set any extra parameters functional
                print('Adding DFT functionals...')
                self.add_dft_functional(singlet_path, functional_name)

                # set symmetry
                print('Updating symmetry...')
                self.set_symmetry(singlet_path, self.symmetry)

                # increase iterlimit
                print('Increasing iterlimit...')
                self.increase_scfiterlimit(singlet_path)

                # copy to triplet, cation, tddft and set occupations
                triplet_path = os.path.join(molname, functional_name, 'triplet')
                cation_path = os.path.join(molname, functional_name, 'cation')
                tddft_path = os.path.join(molname, functional_name, 'tddft')

                if 'triplet' in self.occupation_setup_list:
                    print('Creating triplet...')
                    shutil.copytree(singlet_path, triplet_path)
                    if self.symmetry == 'cs':
                        pass
                    else:
                        self.alle_swap_to_triplet(triplet_path)

                if 'cation' in self.occupation_setup_list:
                    print('Creating cation...')
                    shutil.copytree(singlet_path, cation_path)
                    if self.symmetry == 'cs':
                        pass
                    else:
                        self.alle_swap_to_cation(cation_path)

                if 'tddft' in self.occupation_setup_list:
                    print('Creating TDDFT calc...')
                    shutil.copytree(singlet_path, tddft_path)
                    self.add_tddft_excitation(tddft_path, self.ex_symmetry, self.scfinstab)

                # Now the singlet occ
                self.add_occupation(singlet_path, molecule['total_es'], 'singlet_uhf', self.occ_orb_sym)

                # run calcs if desired
                if self.run_jobs is True:
                    subprocess.call(self.qsub_command, shell=True, cwd=singlet_path)
                    if 'triplet' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=triplet_path)
                    if 'cation' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=cation_path)
                    if 'tddft' in self.occupation_setup_list:
                        subprocess.call(self.qsub_command, shell=True, cwd=tddft_path)

        os.chdir('..')
        return

    def alle_singlet_uhf_setup(self, calc_path):
        print("Setting up all-e UHF singlet...")
        define_cmds_path = 'setup_cmds'

        with open(os.path.join(calc_path, define_cmds_path), 'w') as var_file:
            var_file.writelines(setup_cmds)
        var_file.close()

        command = 'define < %s' % define_cmds_path
        subprocess.call(command, shell=True, cwd=calc_path)
        return

    def alle_swap_to_triplet(self, calc_path):
        print("Changing to triplet...")
        define_cmds_path = 'to_t_cmds'

        with open(os.path.join(calc_path, define_cmds_path), 'w') as var_file:
            var_file.writelines(ref_to_t_cmds)
        var_file.close()

        command = 'define < %s' % define_cmds_path
        subprocess.call(command, shell=True, cwd=calc_path)
        return

    def alle_swap_to_cation(self, calc_path):
        print("Changing to cation...")
        define_cmds_path = 'to_i_cmds'

        with open(os.path.join(calc_path, define_cmds_path), 'w') as var_file:
            var_file.writelines(ref_to_i_cmds)
        var_file.close()

        command = 'define < %s' % define_cmds_path
        subprocess.call(command, shell=True, cwd=calc_path)
        return

    def gather_results(self):
        # gather results. what's the best format?
        # Can output all results into csv. But also can gether them into a list of dicts for each molecule with respect s, t, ie energies.
        # this can be easily manipulated for averages, etc. though I need the all-e input somehow too.
        # Maybe it would be easier just to calculate the all-e mols too?
        # No. Difference of means = mean difference. If I work out the ref average result for a test set, I can subtract the result from my pseudo average

        totals_command = "grep -r \"  total energy\" */*/*/dscf.out > dscf_totals.dat"
        homo_command = "grep -r \"HOMO:\" */*/singlet_uhf/eiger.out > dscf_homos.dat"
        tddft_command = "grep \"Excitation energy:\" */*/*/escf.out > escf_totals.dat"
        subprocess.call(totals_command, shell=True)
        subprocess.call(homo_command, shell=True)
        subprocess.call(tddft_command, shell=True)

        import csv

        #process grep results
        totals_file = open('dscf_totals.dat', 'r')
        totals_file_data = totals_file.readlines()

        homos_file = open('dscf_homos.dat', 'r')
        homos_file_data = homos_file.readlines()

        tddft_file = open('escf_totals.dat', 'r')
        tddft_file_data = tddft_file.readlines()

        data_dict = {}
        csv_rows = [['Molecule', 'Functional', 'Calculation Type', 'Energy (H)']]

        for line in totals_file_data:
            path_part = line.split()[0].split('/')
            molname = path_part[0]
            functional = path_part[1]
            occ_type = path_part[2]
            total_e = line.split()[-2]

            if molname not in data_dict:
                data_dict[molname] = {}

            if functional not in data_dict[molname]:
                data_dict[molname][functional] = {}

            if occ_type not in data_dict[molname][functional]:
                data_dict[molname][functional][occ_type] = total_e

            row = [molname, functional, occ_type, total_e]
            csv_rows.append(row)

        for line in tddft_file_data:
            path_part = line.split()[0].split('/')
            molname = path_part[0]
            functional = path_part[1]
            occ_type = path_part[2]
            ex_e = line.split()[-1]

            row = [molname, functional, occ_type, ex_e]
            data_dict[molname][functional][occ_type] = ex_e
            csv_rows.append(row)

        for line in homos_file_data:
            path_part = line.split()[0].split('/')
            molname = path_part[0]
            functional = path_part[1]
            homo_e = line.split()[-5]

            row = [path_part[0], path_part[1], 'singlet_HOMO', homo_e]
            data_dict[molname][functional]['singlet_HOMO'] = homo_e
            csv_rows.append(row)

        print(data_dict)
        for functional in self.functional_list:
            average_s_e = sum(float(data_dict.get(molecule, 0).get(functional, 0).get('singlet_uhf', 0)) for molecule in data_dict) / len(data_dict)
            average_t_e = sum(float(data_dict.get(molecule, 0).get(functional, 0).get('triplet', 0)) for molecule in data_dict) / len(data_dict)
            average_c_e = sum(float(data_dict.get(molecule, 0).get(functional, 0).get('cation', 0)) for molecule in data_dict) / len(data_dict)
            average_homo_e = sum(float(data_dict.get(molecule, 0).get(functional, 0).get('singlet_HOMO', 0)) for molecule in data_dict) / len(data_dict)
            average_tddft_e = sum(float(data_dict.get(molecule, 0).get(functional, 0).get('tddft', 0)) for molecule in data_dict) / len(data_dict)
            csv_rows.append([functional, 'average singlet E', average_s_e])
            csv_rows.append([functional, 'average triplet E', average_t_e])
            csv_rows.append([functional, 'average cation E', average_c_e])
            csv_rows.append([functional, 'average st gap E', average_t_e - average_s_e])
            csv_rows.append([functional, 'average ionisation E', average_c_e - average_s_e])
            csv_rows.append([functional, 'average HOMO E', average_homo_e])
            csv_rows.append([functional, 'average TDDFT E', average_tddft_e])

        with open('results.csv', 'w') as writeFile:
            writer = csv.writer(writeFile)
            writer.writerows(csv_rows)

        totals_file.close()
        homos_file.close()
        tddft_file.close()
        writeFile.close()

        # parse and insert into csv? Can you copy macros straight from python into csv?
        # needs to end up in latex table format
        # If I can read it sensibly enough into python I could do the averaging and percentage errors manually

        return

if __name__ == "__main__":
    calcrunner = CalculationRunner()
    if 'gather' in sys.argv:
        calcrunner.gather_results()
    elif 'sort' in sys.argv:
        calcrunner.sort_atoms(sys.argv[2])
    elif 'ref' in sys.argv:
        calcrunner.reference_run(sys.argv[2])
    else:
        calcrunner.run(sys.argv[1], sys.argv[2], sys.argv[3])
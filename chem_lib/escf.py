"""This file contains functions to read and handle ESCF output."""
import linecache
import os
import numpy
import subprocess

sample_excitations = [
    {'id': '14a',
     'energy_nm': 150.0,
     'oscillator_vel': 0.0,
     'rotatory_vel': 0.0,
     'dominant_contributions': [
         {
             'occ_orbital': '5a',
             'virt_orbital': '7a',
             'coeff_sq': 46.5,
          }
     ]},
]


class ESCFControl:
    def __init__(self):
        self.calc_folder_path = '.'
        self.escf_file = 'escf.log'
        self.excitation_list = []
        self.dipole_scan_focus = 0.5
        self.undetectable_excitation_replacement_error = 3000.0
        self.undetectable_excitation_oscillator_replacement_error = 1000.0
        self.spyctrum_path = '../spyctrum/spyctrum/spyctrum.py'
        self.convolution_path = 'convolution.csv'

    def read_escf(self):
        """Sorts all excitations into a convenient list with their properties."""
        var_file_path = os.path.join(self.calc_folder_path, '%s.log' % 'escf')
        var_file = open(var_file_path, 'r')
        excitation_dict = {'dominant_contributions': []}

        for i, line in enumerate(var_file):
            split_line = ' '.join(line.split()).split()
            if split_line:
                if split_line[-1] == 'excitation':
                    if 'WARNING!' not in linecache.getline(var_file_path, i + 7):
                        excitation_dict['id'] = split_line[0] + split_line[2]
                        excitation_dict['energy_ev'] = float(' '.join(linecache.getline(var_file_path, i + 8).split()).split()[4])
                        excitation_dict['energy_nm'] = float(' '.join(linecache.getline(var_file_path, i + 10).split()).split()[4])
                        excitation_dict['oscillator_str'] = float(' '.join(linecache.getline(var_file_path, i + 17).split()).split()[2])
                        excitation_dict['rotatory_str'] = float(' '.join(linecache.getline(var_file_path, i + 26).split()).split()[2])
                        linecache.clearcache()
                elif split_line[0] == 'occ.':
                    for dominant_contribution_index in range(1, 8):
                        dominant_contribution_dict = {}
                        split_contribution_line = ' '.join(linecache.getline(var_file_path, i + dominant_contribution_index).split()).split()
                        if len(split_contribution_line) == 7:
                            dominant_contribution_dict['occ_orbital'] = split_contribution_line[0] + split_contribution_line[1]
                            dominant_contribution_dict['virt_orbital'] = split_contribution_line[3] + split_contribution_line[4]
                            dominant_contribution_dict['coeff_sq'] = split_contribution_line[6]
                            excitation_dict['dominant_contributions'].append(dominant_contribution_dict)
                        linecache.clearcache()
                elif 'Electric transition dipole moment (velocity rep.)' in line:
                    excitation_dict['electric_dipole_norm'] = float(
                        ' '.join(linecache.getline(var_file_path, i + 3).split()).split()[3])
                elif 'Electric quadrupole transition moment' in line:
                    self.excitation_list.append(excitation_dict)
                    excitation_dict = {'dominant_contributions': []}

    def run_spyctrum(self, spyctrum_path=self.spyctrum_path, lmin=1, lmax=1200):
        """Runs the spyctrum program to generate the spectral data file."""
        print("Convoluting the problem (spyctrum)...")
        command = '%s -t %s -m convolution -l %s %s' % (self.spyctrum_path,
                                                        self.escf_file,
                                                        lmin,
                                                        lmax)

        subprocess.call(command, cwd=spyctrum_path)

    def read_spyctrum(self, convolution_path):
        var_file_path = os.path.join(convolution_path)
        var_file = open(var_file_path, 'r')
        convolutions = []

        for i, line in enumerate(var_file):
            if line[0] == '#':
                continue
            else:
                convolutions.append(line.split())

        return convolutions

    def evaulate_spyctrum_difference(self):
        """Reads and finds the difference between two spyctrum spectra.
        The spyctra should be comparable i.e. produced over the same wavelength range."""
        reference_convolutions = self.read_spyctrum(self.convolution_path)
        pseudo_convolutions = self.read_spyctrum(self.convolution_path)

        for i, reference_convolution in enumerate(reference_convolutions):
            pseudo_convolutions[i] = pseudo_convolutions[i] - reference_convolution

        total_error = 0.0

        for convolution in pseudo_convolutions:
            total_error += convolution[1]

        return total_error

    def evaulate_spyctrum_hash_difference(self):
        """Possibly stupid idea. Multiply ALL wavelengths and amplitudes from a spyctrum together to get a 'spyctral hash',
        Compare hash. ??? Profit"""
        reference_convolutions = self.read_spyctrum(self.convolution_path)
        pseudo_convolutions = self.read_spyctrum(self.convolution_path)

        for i, reference_convolution in enumerate(reference_convolutions):
            pseudo_convolutions[i] = pseudo_convolutions[i] - reference_convolution

        total_error = 0.0

        for convolution in pseudo_convolutions:
            total_error += convolution[0] * convolution[1]

        return total_error

    def identify_excitation_by_dominant_contributions(self, dominant_contributions):
        for excitation in self.excitation_list:
            try:
                calc_contributions = excitation['dominant_contributions'][0:len(dominant_contributions)]
                if dominant_contributions == [contribution['occ_orbital'] for contribution in calc_contributions]:
                    print('Identified orbital %s as desired orbital' % excitation['id'])
                    return excitation
            except IndexError:
                pass

    def identify_excitation_by_dipole_moment(self, reference_electric_dipole_norm):
        error_list = []
        for excitation in self.excitation_list:
            error_list.append(reference_electric_dipole_norm - excitation['electric_dipole_norm'])

        if error_list:
            if min(error_list) < self.dipole_scan_focus:
                print('Identified excitation %s as a matching excitation' % self.excitation_list[error_list.index(min(error_list))]['id'])
                return self.excitation_list[error_list.index(min(error_list))]

    def identify_excitation_by_id(self, excitation_id):
        for excitation in self.excitation_list:
            if excitation['id'] == excitation_id:
                return excitation

    def evaluate_peak_error(self, reference_wavelength, reference_oscillator_str, reference_electric_dipole_norm, ref_excitation_id=None):
        excitation = self.identify_excitation_by_dipole_moment(reference_electric_dipole_norm)
        #excitation = self.identify_excitation_by_id(ref_excitation_id)
        if excitation:
            wavelength_error = numpy.abs(excitation['energy_ev'] - float(reference_wavelength))
            oscillator_str_error = numpy.abs(excitation['oscillator_str'] - float(reference_oscillator_str))
        else:
            wavelength_error = self.undetectable_excitation_replacement_error
            oscillator_str_error = self.undetectable_excitation_oscillator_replacement_error
        return wavelength_error, oscillator_str_error, excitation

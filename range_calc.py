from chem_lib.optimise import BasisControl


def test_calcs():
    orbital_to_update = {'line_type': '$ecp',
                         'basis_ecp_name': 'h ecp-s',
                         'orbital_descriptor': 's-f'}
    new_functions = [
                    {'exponent': -7.524, 'r^n': 1, 'coefficient': 100.0}
                    ]
    orbital_to_update['functions_list'] = new_functions
    trial_coefficients = [-50.0, -20.0, -10.0, -1.0, -0.1, 0.0, 0.1, 0.2, 0.5, 1.0, 10.0, 100.0, 1000.0, 1200.0, 1500.0]
    trial_coefficients_2 = [250, 500, 750]
    trial_exponents = [100.0, 250.0, 500.0, 550.0, 600.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 750.0, 1000.0, 1500.0, 2000.0]
    trial_coefficients_3 = [x / 1000 for x in range(-3800, -3500, 5)]

    trial_values = trial_exponents

    for coeff in trial_values:
        basis = BasisControl()
        orbital_to_update['functions_list'][0]['coefficient'] = str(coeff)
        basis.update_variable(orbital_to_update)
        basis.run_ridft()

    import linecache

    out_file = open('calcs.txt', 'r')
    output_energies = []
    print 'Checking for energies...'
    for i, line in enumerate(out_file):
        split_line = ' '.join(line.split()).split()
        # if len(split_line) == 6:
        #     if split_line[1] == "total" and split_line[2] == "energy":
        #         total_e = split_line[4]
        #         print trial_values
        #         print "found energy %s" % total_e
        #         output_energies.append(total_e)
        if len(split_line) == 6:
            if split_line[0] == "irrep" and split_line[1] == "1b1u":
                energy_line = ' '.join(linecache.getline('calcs.txt', i+2).split()).split()
                total_e = energy_line[2]
                print trial_values
                print "found energy %s" % total_e
                output_energies.append(total_e)
    print trial_values
    print "energies: %s" % output_energies

    with open('trial_outcomes.txt', 'w') as var_file:
        for index, trial_value in enumerate(trial_values):
            var_file.writelines('%s %s \n' % (trial_value, output_energies[index]))

test_calcs()
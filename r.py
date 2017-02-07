import numpy
import copy
from chem_lib.optimise import BasisControl
from chem_lib.molecular_orbitals import MOSControl


def apply_overlap_integral(a, b, orbital_type):
    """Apply integral result to exponents."""
    p = a+b
    if orbital_type == 's':
        return numpy.pi**(3.0/2.0)/(p**(3.0/2.0))  # s orbitals
    elif orbital_type == 'p':
        return numpy.pi**(3.0/2.0)/(2.0*p**(5.0/2.0))  # p orbitals


def apply_operator_integral(a, b, orbital_type):
    p = a+b
    if orbital_type == 's':
        return 4.0*numpy.pi/(2.0*p**2.0)  # s orbitals
    elif orbital_type == 'p':
        return 4*numpy.pi/(3.0*p**3.0)  # p orbitals


def compute_primary_matrix(alphas, orbital_type, operator=None, betas=None):
    """Generate normalised overlap matrix."""
    if betas == None:
        betas = alphas
    # dim = len(alphas)
    S = numpy.zeros([len(alphas), len(betas)])
    for i in range(len(alphas)):
        for j in range(len(betas)):
            a = alphas[i]
            b = betas[j]
            if operator == 'r':
                S[i, j] = apply_operator_integral(a, b, orbital_type) /\
                          numpy.sqrt(apply_overlap_integral(a, a, orbital_type)
                                     * apply_overlap_integral(b, b, orbital_type))
            else:
                S[i, j] = apply_overlap_integral(a, b, orbital_type) /\
                          numpy.sqrt(apply_overlap_integral(a, a, orbital_type)
                                     * apply_overlap_integral(b, b, orbital_type))
    return S


def contract_matrix(coeff_list, matrix):
    # decide the contraction ratios needed.
    coefficient_matrix = numpy.zeros([len(matrix), len(coeff_list)])

    # create coefficient matrix
    count = 0
    for i, basis_function in enumerate(coeff_list):
        for j, coeff in enumerate(basis_function):
            coefficient_matrix[count][i] = basis_function[j]
            count += 1

    # do matrix multiplication
    # print matrix
    bc = numpy.dot(matrix, coefficient_matrix)
    # print "here", bc
    # print coefficient_matrix
    ab = numpy.dot(coefficient_matrix.transpose(), bc)

    return ab


def normalise_contracted_matrix(contracted_overlap_matrix, contracted_operator_matrix=None):
    if contracted_operator_matrix is None:
        contracted_operator_matrix = contracted_overlap_matrix

    normalised_matrix = numpy.zeros([len(contracted_overlap_matrix), len(contracted_overlap_matrix)])
    for i in range(len(contracted_overlap_matrix)):
        for j in range(len(contracted_overlap_matrix)):
            normalised_matrix[i, j] = contracted_operator_matrix[i, j] / numpy.sqrt(contracted_overlap_matrix[i, i]
                                                                                    * contracted_overlap_matrix[j, j])

    return normalised_matrix


def retrieve_turbomole_information(basis_information, mos_information, basis_file_path='basis', mos_file_path='alpha'):
    # obtain information from basis file
    basis_file_extract = BasisControl()
    basis_file_extract.variable_file_path = basis_file_path
    basis_file_extract.read_variables()
    coeff, alphas = basis_file_extract.get_coefficients_exponents(basis_file_extract.basis_file,
                                                                  basis_information['basis_name'],
                                                                  basis_information['orbital_to_extract'])
    contracted_coefficients, contracted_exponents = basis_file_extract.extract_contracted_coefficients(
        basis_information['basis_name'],
        basis_information['orbital_to_extract'])

    # obtain  information from mos file
    moscontrol = MOSControl()
    moscontrol.variable_file_path = mos_file_path
    moscontrol.read_variables()
    cmos = moscontrol.molecular_orbital_file[mos_information['orbital']][mos_information['coeff_range'][0]:
                                                                         mos_information['coeff_range'][1]]

    print 'Extracted info for %s functions, basis %s with orbital %s and basis functions %s.' % (
        basis_information['orbital_to_extract'],
        basis_information['basis_name'],
        mos_information['orbital'],
        range(mos_information['coeff_range'][0] + 1, mos_information['coeff_range'][1] + 1),
    )

    return {'exponents': alphas,
            'coefficients': coeff,
            'contracted_coefficients': contracted_coefficients,
            'contracted_exponents': contracted_exponents,
            'cmos': cmos,
            'orbital_type': basis_information['orbital_to_extract']}


def calculate_expected_r(tm_info):

    # overlap and expectation matrices
    S_prim = compute_primary_matrix(tm_info['exponents'], tm_info['orbital_type'])
    r_prim = compute_primary_matrix(tm_info['exponents'], tm_info['orbital_type'], operator='r')

    # contract and normalise
    contracted_overlap = contract_matrix(tm_info['contracted_coefficients'], S_prim)
    contracted_operator = contract_matrix(tm_info['contracted_coefficients'], r_prim)
    normalised_operator = normalise_contracted_matrix(contracted_overlap,
                                                      contracted_operator_matrix=contracted_operator)

    # Create vector of MO coefficients to contract the operator matrix
    cmos = numpy.array(tm_info['cmos'])

    bc = numpy.dot(normalised_operator, cmos)
    expectation_value = numpy.dot(cmos.transpose(), bc)

    return expectation_value


def calculate_mo_overlap(tm_info, pp_exponent):
    """Calculates the overlap between given MO and a given PP."""
    contracted_coefficients = tm_info['contracted_coefficients']
    exponents = tm_info['exponents']
    cmos = tm_info['cmos']

    # Calculate overlap of un-normalised primitives and PP.
    overlap_vector = compute_primary_matrix(exponents, tm_info['orbital_type'], betas=[pp_exponent]).flatten()

    # Multiply relevant terms by contraction coefficients
    count = 0
    contracted_vector = copy.deepcopy(contracted_coefficients)
    for i, primitive_coefficients_nest in enumerate(contracted_vector):
        for j in range(len(primitive_coefficients_nest)):
            contracted_vector[i][j] = contracted_vector[i][j] * overlap_vector[count]
            count += 1

    # Multiply relevant terms by MO coefficients
    new_list = []
    for thing in contracted_vector:
        new_list.append(sum(thing))
    mo_overlap = numpy.dot(new_list, cmos)

    return mo_overlap


def main():

    basis_information = {
        'basis_name': 'c def2-SV(P)',
        'orbital_to_extract': 'p'
    }

    mos_ref_information = {
        'orbital': '1 b2',
        'coeff_range': [0, 2]
    }

    mos_calc_information = {
        'orbital': '2 a2"',
        'coeff_range': [0, 2]
    }

    basis_h_information = {
        'basis_name': 'h def-TZVP',
        'orbital_to_extract': 's'
    }

    mos_h_information = {
        'orbital': '1 a',
        'coeff_range': [0, 3]
    }

    basis_c_information = {
        'basis_name': 'c def-SV(P)',
        'orbital_to_extract': 'p'
    }
    mos_c_information = {
        'orbital': '1 t1u',
        'coeff_range': [0, 2]
    }

    basis_h_p_information = {
        'basis_name': 'h def-TZVP',
        'orbital_to_extract': 'p'
    }

    mos_h_p_information = {
        'orbital': '1 t1u',
        'coeff_range': [0, 1]
    }

    tm_info_CH3 = retrieve_turbomole_information(basis_information, mos_ref_information)
    tm_info_H = retrieve_turbomole_information(basis_h_information, mos_h_information, basis_file_path='basis_H', mos_file_path='alpha_H')
    tm_info_Hp = retrieve_turbomole_information(basis_h_p_information, mos_h_p_information, basis_file_path='basis_H_p', mos_file_path='mos_H_p')

    exp_r = calculate_expected_r(tm_info_CH3)

    pp_exponent = 0.295

    mo_overlap = calculate_mo_overlap(tm_info_CH3, pp_exponent)
    print '<r> =', exp_r
    print 'effective total Z =', 5.0/exp_r
    print 'MO PP overlap =', mo_overlap
    print 'dZ from PP =', (5.0/exp_r-1)/mo_overlap**2

if __name__ == "__main__":
    main()

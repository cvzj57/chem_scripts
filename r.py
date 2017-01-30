import numpy
from chem_lib.optimise import BasisControl
from chem_lib.molecular_orbitals import MOSControl


def apply_overlap_integral(a, b):
    """Apply integral result to exponents."""
    p = a+b
    # return numpy.pi**(3.0/2.0)/(p**(3.0/2.0))  # s orbitals
    return numpy.pi**(3.0/2.0)/(2.0*p**(5.0/2.0))  # p orbitals


def apply_operator_integral(a, b):
    p = a+b
    # return 4.0*numpy.pi/(2.0*p**2.0)  # s orbitals
    return 4*numpy.pi/(3.0*p**3.0)  # p orbitals


def compute_primary_matrix(alphas, operator=None, overlap=None):
    """Generate normalised overlap matrix."""
    dim = len(alphas)
    S = numpy.zeros([dim, dim])
    for i in range(dim):
        for j in range(dim):
            a = alphas[i]
            b = alphas[j]
            if operator == 'r':
                S[i, j] = apply_operator_integral(a, b)/numpy.sqrt(apply_overlap_integral(a, a)
                                                                   * apply_overlap_integral(b, b))
            else:
                S[i, j] = apply_overlap_integral(a, b)/numpy.sqrt(apply_overlap_integral(a, a)
                                                                  * apply_overlap_integral(b, b))
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
    bc = numpy.dot(matrix, coefficient_matrix)
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


def calculate_expected_r(basis_information, mos_information, basis_file_path='basis', mos_file_path='alpha'):

    # obtain information from basis file
    basis_file_extract = BasisControl()
    basis_file_extract.variable_file_path = basis_file_path
    basis_file_extract.read_variables()
    coeff, alpha = basis_file_extract.get_coefficients_exponents(basis_file_extract.basis_file,
                                                                 basis_information['basis_name'],
                                                                 basis_information['orbital_to_extract'])
    contracted_coefficients = basis_file_extract.extract_contracted_coefficients(
        basis_information['basis_name'],
        basis_information['orbital_to_extract'])

    # obtain  information from mos file
    moscontrol = MOSControl()
    moscontrol.variable_file_path = mos_file_path
    moscontrol.read_variables()
    cmos = moscontrol.molecular_orbital_file[mos_information['orbital']][mos_information['coeff_range'][0]:
    mos_information['coeff_range'][1]]

    # overlap and expectation matrices
    S_prim = compute_primary_matrix(alpha)
    r_prim = compute_primary_matrix(alpha, operator='r')

    # contract and normalise
    contracted_overlap = contract_matrix(contracted_coefficients, S_prim)
    contracted_operator = contract_matrix(contracted_coefficients, r_prim)
    normalised_operator = normalise_contracted_matrix(contracted_overlap,
                                                      contracted_operator_matrix=contracted_operator)

    # Create vector of MO coefficients to contract the operator matrix
    cmos = numpy.array(cmos)

    bc = numpy.dot(normalised_operator, cmos)
    expectation_value = numpy.dot(cmos.transpose(), bc)

    print 'Calculated <r> for %s functions, basis %s with orbital %s and basis functions %s: %s' % (
        basis_information['orbital_to_extract'],
        basis_information['basis_name'],
        mos_information['orbital'],
        range(mos_information['coeff_range'][0]+1, mos_information['coeff_range'][1]+1),
        expectation_value
    )

    return expectation_value


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

    exp_r = calculate_expected_r(basis_information, mos_ref_information)
    # exp_r = calculate_expected_r(basis_c_information, mos_c_information, basis_file_path='basis_C', mos_file_path='mos_C')
    print 'Z =', 5.0/exp_r
    # calculate_expected_r(basis_information, mos_calc_information, mos_file_path='alpha_2')
    # calculate_expected_r(basis_h_information, mos_h_information, basis_file_path='basis_H', mos_file_path='alpha_H')
    # calculate_expected_r(basis_h_p_information, mos_h_p_information, basis_file_path='basis_H_p', mos_file_path='mos_H_p')


if __name__ == "__main__":
    main()

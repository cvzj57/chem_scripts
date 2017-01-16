import numpy
from chem_lib.optimise import BasisControl
from chem_lib.molecular_orbitals import MOSControl


def apply_overlap_integral(a, b):
    """Apply integral result to exponents."""
    p = a+b
    return numpy.pi**(3.0/2.0)/(8*numpy.sqrt(2.0)*p**(5.0/2.0))


def apply_operator_integral(a, b):
    p = a+b
    return numpy.sqrt(2.0)/(6*numpy.sqrt(numpy.pi)*numpy.sqrt(p))


def get_C_A(basis_file_extract, basis_set):
    """Retrieve coefficients and exponents from basis file."""
    alpha_list=[]
    coeff_list=[]

    def extract_for_shell(basis_function_name):
        sub_alpha_list = []
        sub_coeff_list = []
        for shell in basis_file_extract["$basis"][basis_set].keys():
            if basis_function_name in shell:
                a = basis_file_extract["$basis"][basis_set]
                for exp in a[shell]:
                    sub_alpha_list.append(exp["exponent"])
                    sub_coeff_list.append(exp["coefficient"])
        return sub_alpha_list, sub_coeff_list

    p_alpha, p_coeff = extract_for_shell('p')
    coeff_list.append(p_coeff)
    alpha_list.append(p_alpha)

    coeff_list = [item for sublist in coeff_list for item in sublist]
    alpha_list = [item for sublist in alpha_list for item in sublist]

    return coeff_list, alpha_list


def compute_primary_matrix(alpha, operator=None):
    """Generate normalised overlap matrix."""
    dim=len(alpha)
    S = numpy.zeros([dim,dim])
    for i in range(dim):
        for j in range(dim):
            a=alpha[i]
            b=alpha[j]
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


def calculate_expected_r(basis_information, mos_information):

    # obtain information from basis file
    basis_file_extract = BasisControl()
    coeff, alpha = get_C_A(basis_file_extract.basis_file, basis_information['basis_name'])
    contracted_coefficients = basis_file_extract.extract_coefficients(basis_information['basis_name'],
                                                                      basis_information['orbital_to_extract'])

    # overlap and expectation matrices
    S_prim = compute_primary_matrix(alpha)
    r_prim = compute_primary_matrix(alpha, operator='r')

    # contract and normalise
    contracted_overlap = contract_matrix(contracted_coefficients, S_prim)
    contracted_operator = contract_matrix(contracted_coefficients, r_prim)
    normalised_operator = normalise_contracted_matrix(contracted_overlap,
                                                      contracted_operator_matrix=contracted_operator)

    # obtain  information from mos file
    moscontrol = MOSControl()
    moscontrol.read_variables()
    cmos = moscontrol.molecular_orbital_file[mos_information['orbital']][mos_information['coeff_range'][0]:
                                                                         mos_information['coeff_range'][1]]
    flattened_operator = normalised_operator.flatten()

    overall_coeff_number = 0
    for basis_function_id, list_of_contraction_coeffs in enumerate(contracted_coefficients):
        for function_coeff_i in range(len(list_of_contraction_coeffs)):
            flattened_operator[overall_coeff_number] = flattened_operator[overall_coeff_number] \
                                                       * cmos[basis_function_id]
            overall_coeff_number += 1

    expectation_value = 0
    for i in range(len(flattened_operator)):
        for j in range(len(flattened_operator)):
            expectation_value += flattened_operator[i] * flattened_operator[j] * S_prim[i, j]

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

    mos_information = {
        'orbital': '1 b2',
        'coeff_range': [0, 2]
    }

    calculate_expected_r(basis_information, mos_information)
  
if __name__ == "__main__":
    main()
    
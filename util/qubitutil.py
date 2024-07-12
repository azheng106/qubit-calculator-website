import itertools
from enum import Enum

from sympy import isprime
import numpy as np

from threequbitstate import ThreeQubitState, EqClass


class EntanglementStatus(Enum):
    ENTANGLED = 'State is entangled'
    UNKNOWN = 'State is unknown'
    SEPARABLE = 'State is separable'


def convert_to_basic_state(num_qubits, constant_number) -> str:
    """
    Convert a constant number to its basic state (binary) representation. Ex. C3 => |011>
    """
    num_qubits = int(num_qubits)  # Convert string to int
    if constant_number > (2 ** num_qubits - 1):
        raise ValueError('Constant number >= 2^n')
    binary = bin(constant_number).replace('0b', '')
    while len(binary) < num_qubits:
        binary = '0' + binary
    return binary


def remove_qubit_n(n, pairs_dictionary) -> dict:
    """
    Remove qubit n from given state
    :param n: Qubit number to remove
    :param pairs_dictionary: Constant + basic state pairs that make up the state
    :return: Updated dictionary with qubit n removed
    """
    return {(key[:n] + key[n + 1:]): value for key, value in pairs_dictionary.items()}


def is_entangled_2qubit(pairs_dictionary) -> Enum:
    """
    Check if a 2 qubit state is entangled or separable
    """
    if len(pairs_dictionary) > 4:
        raise ValueError('Not a 2 qubit state')
    c0 = pairs_dictionary.get('00', 0)  # Default to 0 if not found
    c1 = pairs_dictionary.get('01', 0)
    c2 = pairs_dictionary.get('10', 0)
    c3 = pairs_dictionary.get('11', 0)
    if c0 * c3 != c1 * c2:
        return EntanglementStatus.ENTANGLED
    else:
        return EntanglementStatus.SEPARABLE


def is_entangled_3qubit(pairs_dictionary) -> Enum:
    """
    Check if a 3 qubit state is entangled or separable
    """
    if len(pairs_dictionary) > 8:
        raise ValueError('Not a 3 qubit state')
    c0 = pairs_dictionary.get('000', 0)  # Default to 0 if not found
    c1 = pairs_dictionary.get('001', 0)
    c2 = pairs_dictionary.get('010', 0)
    c3 = pairs_dictionary.get('011', 0)
    c4 = pairs_dictionary.get('100', 0)
    c5 = pairs_dictionary.get('101', 0)
    c6 = pairs_dictionary.get('110', 0)
    c7 = pairs_dictionary.get('111', 0)

    state = ThreeQubitState(c0, c1, c2, c3, c4, c5, c6, c7)
    slocc_class = state.get_state_class()
    if slocc_class == EqClass.GHZ or slocc_class == EqClass.W:
        return EntanglementStatus.ENTANGLED
    else:
        return EntanglementStatus.SEPARABLE


def is_entangled(pairs_dictionary, n) -> Enum:
    """
    Check if an n-qubit state is entangled, or separable
    :param pairs_dictionary: State to check
    :param n: Number of qubits in the state
    """
    print(coeff_matrix(pairs_dictionary, n))
    if n == 2:
        return is_entangled_2qubit(pairs_dictionary)

    if n == 3:
        return is_entangled_3qubit(pairs_dictionary)

    if len(pairs_dictionary) == 1:
        return EntanglementStatus.SEPARABLE

    if len(pairs_dictionary) == 4:  # Section 3, m = 4 and n >= 2
        values = list(pairs_dictionary.values())
        for i in range(0, 4):
            for j in range(i+1, 4):
                if i != j and values[i] + values[j] == 1:
                    unchecked_values = [values[k] for k in range(4) if k != i and k != j]
                    if sum(unchecked_values) == 1 and values[i] * values[j] == unchecked_values[0] * unchecked_values[1]:
                        return EntanglementStatus.SEPARABLE

    basic_states = list(pairs_dictionary.keys())
    for i in range(len(basic_states[0])):  # if state has |0> or |1> that can be factored out
        if all(basic_state[i] == '0' for basic_state in basic_states) or \
                all(basic_state[i] == '1' for basic_state in basic_states):
            return EntanglementStatus.SEPARABLE

    if isprime(len(basic_states)):  # m (# non-zero coefficients) is prime
        return EntanglementStatus.ENTANGLED
    else:
        basis_matrix = make_basis_matrix(pairs_dictionary)

        if check_canonical_form(basis_matrix) and check_proportional_rows_and_columns(coeff_matrix(pairs_dictionary, n)):
            return EntanglementStatus.SEPARABLE

    return EntanglementStatus.ENTANGLED


def dict_to_hashable(dictionary):
    """
    Convert a dictionary to a hashable type (tuple) so it can be used as a key in a dictionary
    """
    return tuple(dictionary.items())


def make_basis_matrix(pairs_dictionary) -> np.array:
    """
    Create basis matrix B, where each row of the matrix is a basic state
    """
    keys = list(pairs_dictionary.keys())
    matrix = np.array([list(map(int, key)) for key in keys])
    return matrix


def check_canonical_form(basis_matrix) -> bool:
    """
    Check if a basis matrix (matrix of basis states) can has a permutation that can be converted into canonical form
    """
    for perm in itertools.permutations(range(basis_matrix.shape[1])):
        permuted_matrix = basis_matrix[:, perm]
        for row_perm in itertools.permutations(range(basis_matrix.shape[0])):
            permuted_matrix = permuted_matrix[list(row_perm), :]
            if is_canonical(permuted_matrix):
                return True
    return False


def is_canonical(matrix) -> bool:
    """
    Check if a matrix can be converted into a canonical form (PI1 delta
                                                              PI2 delta)
    """
    rows, cols = matrix.shape
    # Transpose the matrix to work with columns as if they are rows
    matrix = matrix.T
    unique_rows = {tuple(row) for row in matrix}

    if len(unique_rows) != rows:
        return False

    for row in unique_rows:
        sub_matrix = matrix[np.all(matrix == row, axis=1)]
        if not all(np.array_equal(sub_matrix[0], r) for r in sub_matrix):
            return False
    return True


def are_proportional(vec1, vec2):
    """
    Check if two vectors are proportional
    """
    ratio = None
    for i in range(len(vec1)):
        if vec1[i] == 0 and vec2[i] == 0:
            continue
        if vec1[i] == 0 or vec2[i] == 0:
            return False
        if ratio is None:
            ratio = vec1[i] / vec2[i]
        elif vec1[i] / vec2[i] != ratio:
            return False
    return True


def check_proportional_rows_and_columns(coeff_matrix):
    """
    Check if coefficient matrix has proportional rows and columns. Appendix A of notes
    """
    rows, cols = coeff_matrix.shape

    # Check rows for proportionality
    for i in range(rows):
        for j in range(i + 1, rows):
            if not are_proportional(coeff_matrix[i, :], coeff_matrix[j, :]):
                return False
    return True


def coeff_matrix(pairs_dictionary, num_qubits):
    """
    Generate the coefficient matrix from a pairs dictionary. Appendix A of notes
    """
    # Generate all possible basic states for the given number of qubits
    basic_states = [bin(i)[2:].zfill(num_qubits) for i in range(2 ** num_qubits)]

    # Extract coefficients for each state, assume 0 if not present in the dictionary
    coefficients = [pairs_dictionary.get(key, 0) for key in basic_states]

    # Determine the dimensions of the matrix
    dim = int(np.sqrt(len(coefficients)))

    # Reshape the coefficients into a matrix
    matrix = np.array(coefficients).reshape(dim, dim)

    return matrix

print(check_canonical_form(np.array([[0, 0, 0], [0, 1, 1], [1, 0, 0], [1, 1, 1]])))

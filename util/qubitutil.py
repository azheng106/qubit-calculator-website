import itertools
from enum import Enum

from sympy import isprime
import numpy as np

from util.threequbitstate import ThreeQubitState, EqClass


class EntanglementStatus(Enum):
    ENTANGLED = 'State is entangled'
    UNKNOWN = 'State is unknown'
    SEPARABLE = 'State is separable'


def convert_to_basic_state(num_qubits: int, constant_number: int) -> str:
    """
    Convert a constant number to its basic state (binary) representation. Ex. C3 => |011>
    """
    if constant_number > (2 ** num_qubits - 1):
        raise ValueError('Constant number >= 2^n')
    binary = bin(constant_number).replace('0b', '')
    while len(binary) < num_qubits:
        binary = '0' + binary
    return binary


def remove_qubit_n(n: int, pairs_dictionary: dict) -> dict:
    """
    Remove qubit n from given state
    :param n: Qubit number to remove
    :param pairs_dictionary: Constant + basic state pairs that make up the state
    :return: Updated dictionary with qubit n removed
    """
    return {(key[:n] + key[n + 1:]): value for key, value in pairs_dictionary.items()}


def is_entangled_2qubit(pairs_dictionary: dict) -> Enum:
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


def is_entangled_3qubit(pairs_dictionary: dict[str, int]) -> Enum:
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


def is_entangled(pairs_dictionary: dict[str, int], n: int) -> Enum:
    """
    Check if an n-qubit state is entangled, or separable
    :param pairs_dictionary: State to check
    :param n: Number of qubits in the state
    """
    if n == 2:
        return is_entangled_2qubit(pairs_dictionary)

    if n == 3:
        return is_entangled_3qubit(pairs_dictionary)

    if len(pairs_dictionary) == 1:
        return EntanglementStatus.SEPARABLE

    basic_states: list[str] = list(pairs_dictionary.keys())
    if can_be_factored(pairs_dictionary):
        return EntanglementStatus.SEPARABLE

    if isprime(len(basic_states)):  # m (# non-zero coefficients) is prime
        return EntanglementStatus.ENTANGLED

    if len(pairs_dictionary) == 4:  # Section 3, m = 4 and n >= 2
        values = list(pairs_dictionary.values())
        for i in range(0, 4):
            for j in range(i+1, 4):
                if i != j and values[i] + values[j] == 1:
                    unchecked_values = [values[k] for k in range(4) if k != i and k != j]
                    if sum(unchecked_values) == 1 and values[i] * values[j] == unchecked_values[0] * unchecked_values[1]:
                        return EntanglementStatus.SEPARABLE

    basis_matrix = make_basis_matrix(pairs_dictionary)
    if check_canonical_form(basis_matrix) and \
            check_proportional_rows_and_columns(coeff_matrix(pairs_dictionary, n)):
        return EntanglementStatus.SEPARABLE

    return EntanglementStatus.ENTANGLED


def count_num_0(basic_states: list[str], index: int):
    """
    Returns how many 0s are at the index qubit
    Ex. count_num_0(['00', '11', '01', '10'], 0) = 2
    Ex. count_num_0(['00', '11', '10'], 0) = 1
    """
    count = 0
    for state in basic_states:
        if state[index] == '0':
            count += 1
    return count


def can_be_factored(pairs_dict: dict[str, int]) -> bool:
    basic_states = list(pairs_dict.keys())
    num_qubits = len(basic_states[0])
    num_states = len(basic_states)
    check_num_0_count = True
    for i in range(num_qubits):  # if state has |0> or |1> that can be factored out
        if all(basic_state[i] == '0' for basic_state in basic_states) or \
                all(basic_state[i] == '1' for basic_state in basic_states):
            return True
        else:  # Check if can factor in a different way, such as
            # |000> + |001> + |100> + |101> = (|0> + |1>)(|00> + |11>)
            if num_states % 2 != 0:
                continue
            coeffs = list(pairs_dict.values())
            first_half = coeffs[: len(coeffs) // 2]
            second_half = coeffs[len(coeffs) // 2:]
            if first_half != second_half and [n * -1 for n in first_half] != second_half:
                continue
            if check_num_0_count and count_num_0(basic_states, i) == num_states / 2:
                check_num_0_count = False
            else:
                remaining_qubits = set([state[i+1:] for state in basic_states])
                num_unique = len(remaining_qubits)
                if num_unique == num_states / 2:  # and check signs
                    return True
    return False


def make_basis_matrix(pairs_dictionary: dict[str, int]) -> np.array:
    """
    Create basis matrix B, where each row of the matrix is a basic state
    """
    keys = list(pairs_dictionary.keys())
    matrix = np.array([list(map(int, key)) for key in keys])
    return matrix


def check_canonical_form(basis_matrix: np.array) -> bool:
    """
    Check if a basis matrix (matrix of basis states) has a permutation that can be converted into canonical form
    """
    for perm in itertools.permutations(range(basis_matrix.shape[1])):  # all permutations of columns
        permuted_matrix = basis_matrix[:, perm]
        if is_canonical(permuted_matrix):
            return True
    return False


def is_equal_rows(matrix: np.array) -> bool:
    """Check if all rows in the matrix are equal."""
    return np.all(np.all(matrix == matrix[0, :], axis=0))


def is_canonical(matrix: np.array) -> bool:
    """Check if the matrix can be written in canonical form (in any permutation of rows)."""
    num_rows = matrix.shape[0]
    num_cols = matrix.shape[1]
    factorizations = []  # array of tuples

    # Get all factorizations of num_rows
    for i in range(2, num_rows // 2):
        if num_rows % i == 0:
            factorizations.append((i, num_rows // i))

    # Number of PIs and deltas in the canonical form is equal to the first number in the factorization
    # e.g When checking 2x3, there are 2 PIs/deltas. When checking 3x2, there are 3 PIs/deltas
    for factorization in factorizations:
        for i in range(2):  # Check each factorization in both orders, such as 2x3 and 3x2
            deltas = []
            num_pis = factorization[i]
            for barrier_col in range(1, num_cols):  # Define the left-right barrier between PIs and deltas
                left_side = matrix[:, : barrier_col]
                if not is_possible_to_canonical(left_side, num_pis):
                    break  # Saves lots of needless row permutation calculations
                for row_perm in itertools.permutations(range(matrix.shape[0])):
                    matrix = matrix[list(row_perm), :]
                    for pi_height in range(num_rows // num_pis, num_rows + 1, num_rows // num_pis):
                        pi = matrix[pi_height - num_rows // num_pis: pi_height, : barrier_col]
                        if not is_equal_rows(pi):
                            return False
                        delta = matrix[pi_height - num_rows // num_pis: pi_height, barrier_col:]
                        deltas.append(delta)
            # Check if every matrix in deltas is equal to each other
            if not all(np.array_equal(deltas[0], delta) for delta in deltas):
                return False
    return True


def is_possible_to_canonical(matrix: np.array, num_pis: int) -> bool:
    """
    Checks if it's possible for any row permutation of a matrix to be written in canonical form
    It is impossible if
        1. there are more unique rows than the number of pis
        or
        2. there is the incorrect number of each unique row for all rows of each pi to be equal
    :param matrix: The part of basis matrix B that was moved to the left
    :param num_pis: Number of sections that matrix must be split into, with all rows of each individual section equal
    :return: Whether it is possible for any row permutation of matrix to be written in canonical form
    """
    m = matrix.shape[0]
    divisor = m // num_pis
    row_count = {}
    for row in matrix:
        if tuple(row) in row_count:
            row_count[tuple(row)] += 1
        else:
            row_count[tuple(row)] = 1

    for row_tuple, count in row_count.items():
        if count % divisor != 0:
            return False
    if len(row_count.keys()) > num_pis:
        return False
    return True


def are_proportional(vec1: list[int], vec2: list[int]) -> bool:
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


def check_proportional_rows_and_columns(coeff_matrix: np.array) -> bool:
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


def coeff_matrix(pairs_dictionary: dict[str, int], num_qubits: int) -> np.array:
    """
    Generate the coefficient matrix from a pairs dictionary. Appendix A of notes
    """
    # Generate all possible basic states for the given number of qubits
    basic_states = [bin(i)[2:].zfill(num_qubits) for i in range(2 ** num_qubits)]

    # Extract coefficients for each state, assume 0 if not present in the dictionary
    coefficients = [pairs_dictionary.get(key) for key in basic_states if pairs_dictionary.get(key) is not None]

    # Determine the dimensions of the matrix
    n = len(coefficients)
    if np.sqrt(n).is_integer():  # If n is a perfect square
        dim = int(np.sqrt(n))
        matrix = np.array(coefficients).reshape(dim, dim)
    else:  # If n is not a perfect square, find the two closest factors
        factors = [(i, n // i) for i in range(1, int(n**0.5) + 1) if n % i == 0]
        dim = factors[-1]
        matrix = np.array(coefficients).reshape(dim[0], dim[1])

    return matrix

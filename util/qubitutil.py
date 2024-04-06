from enum import Enum

from util.threequbitstate import ThreeQubitState, EqClass


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


cached_results = {}  # Use memoization to drastically reduce run-time of is_entangled() function. Before this, it took 5+ mins for 12-qubit state. Now it takes a few seconds for 31 qubits!


def is_entangled(pairs_dictionary, n) -> Enum:
    """
    Check if an n-qubit state is entangled, unknown, or separable. A state is entangled if two of its children (with one qubit removed) are entangled
    :param pairs_dictionary: State to check
    :param n: Number of qubits in the state
    """

    if n == 2:
        return is_entangled_2qubit(pairs_dictionary)

    if n == 3:
        return is_entangled_3qubit(pairs_dictionary)

    if len(pairs_dictionary) == 1:
        return EntanglementStatus.SEPARABLE

    if dict_to_hashable(pairs_dictionary) in cached_results:  # If result has previously been calculated, use cached result instead of calculating it again.
        return cached_results.get(dict_to_hashable(pairs_dictionary))

    basic_states = list(pairs_dictionary.keys())
    for i in range(len(basic_states[0])):
        if all(basic_state[i] == '0' for basic_state in basic_states) or \
                all(basic_state[i] == '1' for basic_state in basic_states):
            cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.SEPARABLE
            return EntanglementStatus.SEPARABLE

    entangled_count = 0
    for i in range(n):
        qubit_i_removed = remove_qubit_n(i, pairs_dictionary)
        if is_entangled(qubit_i_removed, n - 1) == EntanglementStatus.ENTANGLED:
            entangled_count += 1
        if entangled_count >= 2:
            cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.ENTANGLED
            return EntanglementStatus.ENTANGLED

    cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.UNKNOWN
    return EntanglementStatus.UNKNOWN


def dict_to_hashable(dictionary):
    """
    Convert a dictionary to a hashable type (tuple) so it can be used as a key in a dictionary
    """
    return tuple(dictionary.items())

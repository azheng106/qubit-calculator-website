from enum import Enum
from sympy import isprime


class EntanglementStatus(Enum):
    ENTANGLED = 'State is entangled'
    UNKNOWN = 'State is unknown'
    SEPARABLE = 'State is separable'


def convert_to_basic_state(num_qubits, constant_number) -> str:
    """
    Convert a constant number to its basic state (binary) representation. Ex. C3 is |011>
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
    Check if a 2 qubit state is entangled
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


cached_results = {}  # Use memoization to drastically reduce run-time of is_entangled() function. Before this, it took 5+ mins for 12-qubit state. Now it takes a few seconds for 31 qubits!


def is_entangled(pairs_dictionary, n) -> Enum:
    """
    Check if an n-qubit state is entangled, unknown, or separable. A state is entangled if two of its children (with one qubit removed) are entangled
    :param pairs_dictionary: State to check
    :param n: Number of qubits in the state
    """
    # print(f'State: {pairs_dictionary}')

    if n == 2:
        return is_entangled_2qubit(pairs_dictionary)

    if dict_to_hashable(pairs_dictionary) in cached_results:  # If result has previously been calculated, use cached result instead of calculating it again.
        return cached_results.get(dict_to_hashable(pairs_dictionary))

    if isprime(len(pairs_dictionary)):
        basic_states = list(pairs_dictionary.keys())
        for i in range(len(basic_states[0])):
            if all(basic_state[i] == '0' for basic_state in basic_states) or \
                    all(basic_state[i] == '1' for basic_state in basic_states):
                cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.SEPARABLE
                return EntanglementStatus.SEPARABLE
        cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.ENTANGLED
        return EntanglementStatus.ENTANGLED
    else:
        if n == 2:
            return is_entangled_2qubit(pairs_dictionary)
        else:
            entangled_count = 0
            for i in range(n):
                qubit_i_removed = remove_qubit_n(i, pairs_dictionary)
                if is_entangled(qubit_i_removed, n - 1) == EntanglementStatus.ENTANGLED:
                    entangled_count += 1
                    # print(f'{qubit_i_removed} is entangled. Entangled Count for {pairs_dictionary}: {entangled_count}')
                if entangled_count >= 2:
                    # print(f'Entangled count >2, {pairs_dictionary} is entangled')
                    cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.ENTANGLED
                    return EntanglementStatus.ENTANGLED
        cached_results[dict_to_hashable(pairs_dictionary)] = EntanglementStatus.UNKNOWN
        return EntanglementStatus.UNKNOWN


def dict_to_hashable(dictionary):
    """
    Convert a dictionary to a hashable type (tuple) so it can be used as a key in a dictionary
    """
    return tuple(dictionary.items())


""" Old calculation without prime number thing
def old_calc(pairs_dictionary, n) -> Enum:
    if n == 2:
        return is_entangled_2qubit(pairs_dictionary)
    else:
        entangled_count = 0
        for i in range(n):
            qubit_i_removed = remove_qubit_n(i, pairs_dictionary)
            if old_calc(qubit_i_removed, n - 1) == EntanglementStatus.ENTANGLED:
                entangled_count += 1
                # print(f'{qubit_i_removed} is entangled. Entangled Count for {pairs_dictionary}: {entangled_count}')
            if entangled_count >= 2:
                # print(f'Entangled count >2, {pairs_dictionary} is entangled')
                return EntanglementStatus.ENTANGLED
    return EntanglementStatus.UNKNOWN
"""

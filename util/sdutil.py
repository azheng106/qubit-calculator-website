import cmath

import numpy as np

# Schmidt Decomposition useful matrices
identity = np.identity(2)
pauliX = np.array([[0, 1], [1, 0]])
pauliZ = np.array([[1, 0], [0, -1]])


def calculate_UA(alpha, beta, gamma, tau) -> list:
    U_A = []
    if alpha == 0 and beta == 0 and gamma == 0 and tau == 0:  # A-B-C, B-AC, C-AB
        U_A.append(identity)
        print("Case 1")
    elif alpha == 0 and gamma == 0 and beta != 0:  # GHZ
        U_A.append(identity)
        U_A.append(pauliX)
        print("Case 2")
    elif alpha == 0 and beta == 0 and gamma != 0:  # W, A-BC
        U_A.append(identity)
        print("Case 3")
    elif alpha == 0 and beta != 0 and gamma != 0:  # GHZ
        delta = abs(beta) / cmath.sqrt(abs(beta) ** 2 + abs(gamma) ** 2)
        U_A.append(identity)
        U_A.append(np.array([[-gamma * delta / beta, delta], [delta, np.conj(gamma) * delta / np.conj(beta)]]))
        print("Case 4")
    elif alpha != 0 and beta == 0 and gamma == 0:  # W, A-BC
        U_A.append(pauliX)
        print("Case 5")
    elif alpha != 0 and gamma == 0 and beta != 0:  # GHZ
        delta = abs(beta) / cmath.sqrt(beta ** 2 + abs(alpha) ** 2)
        U_A.append(pauliX)
        U_A.append(np.array([[delta, -alpha * delta / beta], [-np.conj(alpha) * delta / np.conj(beta), delta]]))
        print("Case 6")
    elif alpha != 0 and gamma != 0:
        if beta == 0:  # GHZ
            t = cmath.sqrt(-alpha * gamma) / alpha
            U_A.append(big_T_array(t))
            U_A.append(big_T_array(-t))
            print("Case 7.1")
        elif beta != 0 and tau == 0:  # W, A-BC
            t = -beta / (2 * alpha)
            U_A.append(big_T_array(t))
            print("Case 7.2")
        elif beta != 0 and tau != 0:  # GHZ
            t = (-beta - cmath.sqrt(tau)) / (2 * alpha)
            U_A.append(big_T_array(t))
            t = (-beta + cmath.sqrt(tau)) / (2 * alpha)
            U_A.append(big_T_array(t))
            print("Case 7.3")
    return U_A


def big_T_array(t) -> np.array:
    denom = cmath.sqrt(abs(t) ** 2 + 1)
    return np.array([[t / denom, 1 / denom],
                     [1 / denom, -np.conj(t) / denom]])

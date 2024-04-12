import cmath
from enum import Enum

import numpy as np

from util.evalutil import apply_tolerance
from util.sdutil import calculate_UA, identity, pauliX, pauliZ


class ThreeQubitState:
    # c0 * |000> + c1 * |001> + ... + c7 * |111>
    def __init__(self, c0, c1, c2, c3, c4, c5, c6, c7):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        self.c6 = c6
        self.c7 = c7

    def get_3tangle(self) -> float:
        """
        Returns the 3-tangle value of the state, which is used in determining the SLOCC eq. class
        """
        d = ((self.c0 * self.c7 + self.c1 * self.c6 - self.c2 * self.c5 - self.c3 * self.c4) ** 2
             - 4 * (self.c2 * self.c4 - self.c0 * self.c6) * (self.c3 * self.c5 - self.c1 * self.c7))
        tau = 4 * abs(d)
        return tau

    def get_state_class(self) -> Enum:
        """
        Returns the SLOCC equivalence class of the 3-qubit state, based on its tensor decomposition
        """
        eq10 = (self.c0 * self.c7 + self.c1 * self.c6 - self.c2 * self.c5 - self.c3 * self.c4) == 0
        eq11 = (self.c0 * self.c7 - self.c1 * self.c6 - self.c2 * self.c5 + self.c3 * self.c4) == 0

        test1 = (self.c0 * self.c3 == self.c1 * self.c2) and (self.c5 * self.c6 == self.c4 * self.c7)
        test2 = (self.c1 * self.c4 == self.c0 * self.c5) and (self.c3 * self.c6 == self.c2 * self.c7)
        test3 = (self.c3 * self.c5 == self.c1 * self.c7) and (self.c2 * self.c4 == self.c0 * self.c6)

        if eq10:
            if test1 and test2 and test3:
                return EqClass.A_B_C

            if not test1 and test2 and test3:
                return EqClass.A_BC

            if test1 and not test2 and test3:
                return EqClass.B_AC

        if eq11:
            if test1 and test2 and not test3:
                return EqClass.C_AB

        if self.get_3tangle() == 0:
            if not test1 and not test2 and not test3:
                return EqClass.W
        else:
            return EqClass.GHZ

    def get_schmidt_decomposition(self) -> list:
        alpha = apply_tolerance(self.c0 * self.c3 - self.c1 * self.c2)
        beta = apply_tolerance(self.c0 * self.c7 - self.c1 * self.c6 - self.c2 * self.c5 + self.c3 * self.c4)
        gamma = apply_tolerance(self.c4 * self.c7 - self.c5 * self.c6)
        tau = apply_tolerance(self.get_3tangle() / 4)

        C_0 = np.array([[self.c0, self.c1], [self.c2, self.c3]])
        C_1 = np.array([[self.c4, self.c5], [self.c6, self.c7]])

        U_A = calculate_UA(alpha, beta, gamma, tau)
        schmidt_decompositions = []

        if len(U_A) == 0:
            raise Exception('0 returned values for U^A, calculation error')
        for UA in U_A:
            """
            Appendix E, calculate SVD for degenerate (det=0) matrix L_0A
            Calculate C_0_prime, which is the first half of Schmidt Decomposition
            """
            U_B = None
            U_C = None
            C_0_prime = None
            L_0A = apply_tolerance(UA[0, 0] * C_0 + UA[0, 1] * C_1)
            if np.linalg.det(L_0A) != 0:
                raise Exception('L_0A determinant != 0')

            if np.all(L_0A == 0):  # Case 1
                print("1")
                U_B = identity
                U_C = identity
                C_0_prime = np.array([[0, 0], [0, 0]])
            elif np.count_nonzero(L_0A) == 1:  # Case 2
                print("2")
                nonzero_loc = np.transpose(np.nonzero(L_0A))
                row_ind, col_ind = nonzero_loc[0]
                r_exp_i_phi = L_0A[row_ind, col_ind]
                r = abs(r_exp_i_phi)
                exp_i_phi = r_exp_i_phi / r
                U = np.array([[exp_i_phi ** -1, 0], [0, exp_i_phi]])
                U_prime = np.array([[0, exp_i_phi ** -1], [exp_i_phi, 0]])
                sigma = np.array([[r, 0], [0, 0]])
                if tuple(nonzero_loc[0]) == (0, 0):
                    U_B = U
                    U_C = identity
                    C_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (0, 1):
                    print("2.2")
                    U_B = U
                    U_C = pauliX
                    C_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (1, 0):
                    U_B = U_prime
                    U_C = identity
                    C_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (1, 1):
                    U_B = U_prime
                    U_C = pauliX
                    C_0_prime = sigma
            elif np.count_nonzero(L_0A) == 2:  # Case 3
                print("3")
                pos1, pos2 = np.transpose(np.nonzero(L_0A))
                a = L_0A[pos1[0], pos1[1]]
                b = L_0A[pos2[0], pos2[1]]
                temp_a_b = cmath.sqrt(a * np.conj(a) + b * np.conj(b))  # commonly used expression in matrices PI and D
                PI = np.array([[np.conj(a) / temp_a_b, np.conj(b) / temp_a_b],
                               [np.conj(b) / temp_a_b, -(a * np.conj(b) / b) / temp_a_b]])
                D = np.array([[temp_a_b, 0], [0, 0]])
                if np.all(pos1 == (0, 0)) and np.all(pos2 == (0, 1)):  # 3.1
                    U_B = pauliZ
                    U_C = PI
                    C_0_prime = D
                elif np.all(pos1 == (1, 0)) and np.all(pos2 == (1, 1)):  # 3.2
                    U_B = pauliX
                    U_C = PI
                    C_0_prime = D
                elif np.all(pos1 == (0, 0)) and np.all(pos2 == (1, 0)):  # 3.3
                    U_B = PI
                    U_C = identity
                    C_0_prime = D
                elif np.all(pos1 == (0, 1)) and np.all(pos2 == (1, 1)):  # 3.4
                    U_B = PI
                    U_C = pauliX
                    C_0_prime = D
            elif np.count_nonzero(L_0A) == 4:  # Case 4
                print("4")
                a = L_0A[0, 0]
                b = L_0A[0, 1]
                ka = L_0A[1, 0]
                k = ka / a
                temp_k = cmath.sqrt(k * np.conj(k) + 1)
                K = np.array([[1 / temp_k, np.conj(k) / temp_k], [k / temp_k, -1 / temp_k]])

                temp_a_b = cmath.sqrt(a * np.conj(a) + b * np.conj(b))  # Same PI array and temp_a_b as case 3
                PI = np.array([[np.conj(a) / temp_a_b, np.conj(b) / temp_a_b],
                               [np.conj(b) / temp_a_b, (-a * np.conj(b) / b) / temp_a_b]])
                U_B = K
                U_C = PI
                C_0_prime = np.array([[temp_a_b * temp_k, 0], [0, 0]])
            else:  # L_0A cannot have 3 nonzero entries b/c determinant must = 0
                raise Exception('L_0A has 3 nonzero entries. error')

            """
            Calculate C_1_prime, second half of Schmidt Decomposition
            """
            L_1A = UA[1, 0] * C_0 + UA[1, 1] * C_1
            if any(var is None for var in (U_B, U_C, C_0_prime)):
                raise RuntimeError('U_B, U_C, or C_0_prime was not computed')
            C_1_prime = U_B @ L_1A @ np.transpose(U_C)  # @ means matrix multiplication
            sd_form = [C_0_prime[0, 0], C_1_prime[0, 0], C_1_prime[0, 1], C_1_prime[1, 0],
                       C_1_prime[1, 1]]  # Coefficients for |000>, |100>, |101>, |110>, |111>, Eq. 18

            # 1.2
            phases = [cmath.phase(sd_form[i]) for i in range(1, 5)]

            phi = phases[0] - phases[1] - phases[2] + phases[3]
            sd_form = [abs(sd_form[i]) for i in range(0, 5)]
            sd_form[1] *= cmath.exp(1j * phi)  # Multiply |100> coeff by e^(i*phi)

            schmidt_decompositions.append(sd_form)
        if len(schmidt_decompositions) == 2 and schmidt_decompositions[0] == schmidt_decompositions[1]:  # Remove duplicate if both SDs are equal
            schmidt_decompositions.pop()
        return schmidt_decompositions


class EqClass(Enum):
    A_B_C = 'A-B-C'
    A_BC = 'A-BC'
    C_AB = 'C-AB'
    B_AC = 'B-AC'
    W = 'W'
    GHZ = 'GHZ'

import cmath
from enum import Enum

import numpy as np

from util.evalutil import apply_tolerance, truncate, lists_are_equal
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
        tau = apply_tolerance(beta ** 2 - 4 * alpha * gamma)

        c_0 = np.array([[self.c0, self.c1], [self.c2, self.c3]])
        c_1 = np.array([[self.c4, self.c5], [self.c6, self.c7]])

        u_a = calculate_UA(alpha, beta, gamma, tau)
        schmidt_decompositions = []

        if len(u_a) == 0:
            raise Exception('0 returned values for u^A, calculation error')
        for ua in u_a:
            """
            Appendix E, calculate SVD for degenerate (det=0) matrix L_0A
            Calculate c_0_prime, which is the first half of Schmidt Decomposition
            """
            u_b = None
            u_c = None
            c_0_prime = None
            l_0_a = apply_tolerance(ua[0, 0] * c_0 + ua[0, 1] * c_1)

            if np.all(l_0_a == 0):  # Case 1
                u_b = identity
                u_c = identity
                c_0_prime = np.array([[0, 0], [0, 0]])
            elif np.count_nonzero(l_0_a) == 1:  # Case 2
                nonzero_loc = np.transpose(np.nonzero(l_0_a))
                row_ind, col_ind = nonzero_loc[0]
                r_exp_i_phi = l_0_a[row_ind, col_ind]
                r = abs(r_exp_i_phi)
                exp_i_phi = r_exp_i_phi / r
                u = np.array([[exp_i_phi ** -1, 0], [0, exp_i_phi]])
                u_prime = np.array([[0, exp_i_phi ** -1], [exp_i_phi, 0]])
                sigma = np.array([[r, 0], [0, 0]])
                if tuple(nonzero_loc[0]) == (0, 0):
                    u_b = u
                    u_c = identity
                    c_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (0, 1):
                    u_b = u
                    u_c = pauliX
                    c_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (1, 0):
                    u_b = u_prime
                    u_c = identity
                    c_0_prime = sigma
                elif tuple(nonzero_loc[0]) == (1, 1):
                    u_b = u_prime
                    u_c = pauliX
                    c_0_prime = sigma
            elif np.count_nonzero(l_0_a) == 2:  # Case 3
                pos1, pos2 = np.transpose(np.nonzero(l_0_a))
                a = l_0_a[pos1[0], pos1[1]]
                b = l_0_a[pos2[0], pos2[1]]
                temp_a_b = cmath.sqrt(a * np.conj(a) + b * np.conj(b))  # commonly used expression in matrices pi and d
                pi = np.array([[np.conj(a) / temp_a_b, np.conj(b) / temp_a_b],
                               [np.conj(b) / temp_a_b, -(a * np.conj(b) / b) / temp_a_b]])
                d = np.array([[temp_a_b, 0], [0, 0]])
                if np.all(pos1 == (0, 0)) and np.all(pos2 == (0, 1)):  # 3.1
                    u_b = pauliZ
                    u_c = pi
                    c_0_prime = d
                elif np.all(pos1 == (1, 0)) and np.all(pos2 == (1, 1)):  # 3.2
                    u_b = pauliX
                    u_c = pi
                    c_0_prime = d
                elif np.all(pos1 == (0, 0)) and np.all(pos2 == (1, 0)):  # 3.3
                    u_b = pi
                    u_c = identity
                    c_0_prime = d
                elif np.all(pos1 == (0, 1)) and np.all(pos2 == (1, 1)):  # 3.4
                    u_b = pi
                    u_c = pauliX
                    c_0_prime = d
            elif np.count_nonzero(l_0_a) == 4:  # Case 4
                a = l_0_a[0, 0]
                b = l_0_a[0, 1]
                ka = l_0_a[1, 0]
                k = ka / a
                temp_k = cmath.sqrt(k * np.conj(k) + 1)
                K = np.array([[1 / temp_k, np.conj(k) / temp_k], [k / temp_k, -1 / temp_k]])

                temp_a_b = cmath.sqrt(a * np.conj(a) + b * np.conj(b))  # Same pi array and temp_a_b as case 3
                pi = np.array([[np.conj(a) / temp_a_b, np.conj(b) / temp_a_b],
                               [np.conj(b) / temp_a_b, (-a * np.conj(b) / b) / temp_a_b]])
                u_b = K
                u_c = pi
                c_0_prime = np.array([[temp_a_b * temp_k, 0], [0, 0]])
            else:  # L_0A cannot have 3 nonzero entries b/c determinant must = 0
                raise Exception('L_0A has 3 nonzero entries. error')

            """
            Calculate c_1_prime, second half of Schmidt Decomposition
            """
            l_1_a = ua[1, 0] * c_0 + ua[1, 1] * c_1
            if any(var is None for var in (u_b, u_c, c_0_prime)):
                raise RuntimeError('u_b, u_c, or c_0_prime was not computed')
            c_1_prime = u_b @ l_1_a @ np.transpose(u_c)  # @ means matrix multiplication
            sd_form = [c_0_prime[0, 0], c_1_prime[0, 0], c_1_prime[0, 1], c_1_prime[1, 0],
                       c_1_prime[1, 1]]  # Coefficients for |000>, |100>, |101>, |110>, |111>, Eq. 18

            # 1.2
            phases = [cmath.phase(sd_form[i]) for i in range(1, 5)]
            phi = phases[0] - phases[1] - phases[2] + phases[3]
            sd_form = [abs(sd_form[i]) for i in range(0, 5)]
            sd_form[1] *= cmath.exp(1j * phi)  # Multiply |100> coeff by e^(i*phi)

            schmidt_decompositions.append(sd_form)

        if len(schmidt_decompositions) == 2 and lists_are_equal(schmidt_decompositions[0], schmidt_decompositions[1]):  # Remove duplicate if both SDs are equal
            schmidt_decompositions.pop()
        
        schmidt_decompositions = [[truncate(num, 4) for num in inner_list] for inner_list in schmidt_decompositions]  # Truncate trailing decimals

        """
        Convert 2nd term in SD from rectangular (a+bi) to polar (r*e^{ix}) form
        Create an HTML string so website can display e^{ix} with proper superscript notation for exponent
        """
        for i in range(len(schmidt_decompositions)):
            sd = schmidt_decompositions[i]
            r, phi = cmath.polar(sd[1])
            r, phi = truncate(r, 3), truncate(phi, 3)  # Truncate trailing decimals
            if phi == 0:  # So program displays r instead of r * e^(i*0.0)
                continue
            schmidt_decompositions[i] = f'[{sd[0]}, {r}e<sup>{phi}i</sup>, {sd[2]}, {sd[3]}, {sd[4]}]'
        return schmidt_decompositions


class EqClass(Enum):
    A_B_C = 'A-B-C'
    A_BC = 'A-BC'
    C_AB = 'C-AB'
    B_AC = 'B-AC'
    W = 'W'
    GHZ = 'GHZ'

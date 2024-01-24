import cmath
from enum import Enum


class State:
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
        tau = 4 * abs(d)  # If dealing with complex numbers , we'll probably have to change this
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

    def is_normal(self) -> bool:
        x = self.c0**2
        y = self.c7**2
        return (self.c0 ** 2 + self.c1 ** 2 + self.c2 ** 2 + self.c3 ** 2 + self.c4 ** 2 + self.c5 ** 2 + self.c6 ** 2 + self.c7 ** 2) == 1


class EqClass(Enum):
    A_B_C = 'A-B-C'
    A_BC = 'A-BC'
    C_AB = 'C-AB'
    B_AC = 'B-AC'
    W = 'W'
    GHZ = 'GHZ'

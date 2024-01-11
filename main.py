from enum import Enum
import math


class State:
    def __init__(self, c0, c1, c2, c3, c4, c5, c6, c7):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        self.c6 = c6
        self.c7 = c7

    def get_3tangle(self):
        d = ((self.c0 * self.c7 + self.c1 * self.c6 - self.c2 * self.c5 - self.c3 * self.c4) ** 2
             - 4(self.c2 * self.c4 - self.c0 * self.c6)(self.c3 * self.c5 - self.c1 * self.c7))
        tau = 4 * math.abs(d)  # If dealing with complex numbers , we'll probably have to change this
        return tau

    def get_state_class(self):
        # eq10
        eq10 = (self.c0 * self.c7 + self.c1 * self.c6 - self.c2 * self.c5 - self.c3 * self.c4) == 0
        eq11 = (self.c0 * self.c7 - self.c1 * self.c6 - self.c2 * self.c5 + self.c3 * self.c4) == 0

        test1 = (self.c0 * self.c3 == self.c1 * self.c2) and (self.c5 * self.c6 == self.c4 * self.c7)
        test2 = (self.c1 * self.c4 == self.c0 * self.c5) and (self.c3 * self.c6 == self.c2 * self.c7)
        test3 = (self.c3 * self.c5 == self.c1 * self.c7) and (self.c2 * self.c4 == self.c0 * self.c6)
        # A-B-C, A-BC, and B-AC
        if eq10:
            # A-B-C
            if test1 and test2 and test3:
                return EqClasses.A_B_C

            # A-BC
            if not test1 and test2 and test3:
                return EqClasses.A_BC
            # B-AC
            if test1 and not test2 and test3:
                return EqClasses.B_AC

        # C-AB
        if eq11:
            if test1 and test2 and not test3:
                return EqClasses.C_AB

        # For W and GHZ
        # Firstly, W
        if self.get_3tangle() == 0:
            if not test1 and not test2 and not test3:
                # i stg i thought you could use ! instead of not in python but it gives me syntax error
                return EqClasses.W
        # For GHZ
        else:
            return EqClasses.GHZ


# oops i forgot to use your enum but i gotta go so imma just commit


class EqClasses(Enum):
    A_B_C = 1
    A_BC = 2
    C_AB = 3
    B_AC = 4
    W = 5
    GHZ = 6

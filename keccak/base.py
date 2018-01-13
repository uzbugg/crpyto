from abc import ABCMeta, abstractmethod
from math import log10, pow


class Keccak(metaclass=ABCMeta):
    '''
    All parameters based on user input l.
    Every variable with '_' prepending it is private.

    Parameters and algorithms developed according to [1]

    References:
    [1] https://keccak.team/files/Keccak-reference-3.0.pdf
    '''

    # Lane width
    _w = None

    # Permutation width
    _b = None

    # Number of rounds
    _n = None

    # State array
    _a = None

    # Bit array
    _s = None

    # State array dim x
    _x = 5

    # State array dim y
    _y = 5

    # Lane widths
    _lws = [1, 2, 4, 8, 16, 32, 64]

    # Square matrix (Rho/Pi usage)
    _sm = [[0, 1],
           [2, 3]]

    # Column vector (Rho usage)
    _cv = [1, 0]

    # Possible t values
    _t = [i for i in range(24)]

    # Round constants
    # Default insantiation based on _w: 64
    _RC = [0x0000000000000001,
           0x0000000000008082,
           0x800000000000808A,
           0x8000000080008000,
           0x000000000000808B,
           0x0000000080000001,
           0x8000000080008081,
           0x8000000000008009,
           0x000000000000008A,
           0x0000000000000088,
           0x0000000080008009,
           0x000000008000000A,
           0x000000008000808B,
           0x800000000000008B,
           0x8000000000008089,
           0x8000000000008003,
           0x8000000000008002,
           0x8000000000000080,
           0x000000000000800A,
           0x800000008000000A,
           0x8000000080008081,
           0x8000000000008080,
           0x0000000080000001,
           0x8000000080008008
           ]

    def __init__(self, l):
        try:

            # Obtain lane width
            self._w = self._lws[l]

            # Compute perm width
            self._b = 25 * self._w

            # Set state array size
            self._a = [[[0] * self._x] * self._y] * self._w

            # Set bit array size
            self._s = [0] * self._b - 1

            # Get number of rounds
            self._n = 12 + 2 * l

        except IndexError:
            pass

    def keccak_f(self):
        for i in self._n - 1:
            self._a = round(self._a, self.RC[i])
        return self._a

    # Theta auxiliary method
    def sum_a(self, x, z):
        _sum_a = self._a.copy()
        for y in range(4):
            _sum_a += self._a[(x - 1) % 5][y % 5][z % self._w]
        return _sum_a

    # Theta auxiliary method
    def sum_b(self, x, z):
        _sum_b = self._a.copy()
        for y in range(4):
            _sum_b += self._a[(x + 1) % 5][y % 5][(z - 1) % self._w]
        return _sum_b

    # Rho auxiliary method
    def t(self, x, y):
        """
        _sm index map
          idx : val
        -------------
           00 : 0
           01 : 1
           10 : 2
           11 : 3
        """
        tmp = [[None, None],
               [None, None]]

        for t in self._t:

            # Compute index 00 for _tmp, GF(5)
            tmp[0][0] = (pow(((((self._sm[0][0] * self._sm[0][0]) % 5) +
                               (self._sm[0][1] * self._sm[1][0]) % 5) % 5),
                             t) % 5)

            # Compute index 01 for _tmp, GF(5)
            tmp[0][1] = (pow(((((self._sm[0][0] * self._sm[0][1]) % 5) +
                               (self._sm[0][1] * self._sm[1][1]) % 5) % 5),
                             t) % 5)

            # Compute index 10 for _tmp, GF(5)
            tmp[1][0] = (pow(((((self._sm[1][0] * self._sm[0][0]) % 5) +
                               (self._sm[1][1] * self._sm[1][0]) % 5) % 5),
                             t) % 5)

            # Compute index 11 for _tmp, GF(5)
            tmp[1][1] = (pow(((((self._sm[1][0] * self._sm[0][1]) % 5) +
                               (self._sm[1][1] * self._sm[1][1]) % 5) % 5),
                             t) % 5)

            if (((tmp[0][0] * self._cv[0]) % 5 +
                 (tmp[0][1] * self._cv[1]) % 5) % 5 == x and
                ((tmp[1][0] * self._cv[0]) % 5 +
                 (tmp[1][1] * self._cv[1]) % 5) % 5 == y):

                break

        return t
 
    def round(self):

        # Theta
        for x in range(self._x):
            for y in range(self._y):
                for z in range(self._w):
                    self._a[x][y][z] = (self._a[x][y][z] +
                                        (self.sum_a(x, z) +
                                        self.sum_b(x, z)))

        # Rho
        for x in range(self._x):
            for y in range(self._y):
                for z in range(self._w):
                    if x == y == 0:
                        t = -1
                    else:
                        t = self.t()
                    self._a[x][y][z] += self._a[x][y][(z-(t+1)(t+2)/2) % self._w]





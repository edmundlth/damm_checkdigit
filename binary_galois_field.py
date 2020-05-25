import numpy as np

from irreducible_binary_polynomial import IRREDUCIBLE_BINARY_POLYNOMIALS


class BinaryGaloisField(object):
    def __init__(self, degree: int):
        self.degree = degree
        self.order = 2 ** degree
        self.irreducible_polynomial = int(np.sum([2 ** a for a in IRREDUCIBLE_BINARY_POLYNOMIALS[degree]])) + 1

    def quotient_map(self, k):
        result = k
        bit_length = result.bit_length()
        while result >= self.order:
            result ^= (self.irreducible_polynomial << (bit_length - self.degree - 1))
            bit_length = result.bit_length()
        return result

    def mult(self, i, j):
        self._is_elem(i)
        self._is_elem(j)
        result = 0
        bit_string = bin(j)[2:]
        bit_length = j.bit_length()
        for bit_index in range(bit_length):
            bit = bit_string[bit_index]
            result ^= (i * int(bit)) << (bit_length - bit_index - 1)
        return self.quotient_map(result)

    def add(self, i, j):
        self._is_elem(i)
        self._is_elem(j)
        return i ^ j

    def inv(self, k):
        for i in range(self.order):
            if self.mult(i, k) == 1:
                return i
        return

    def div(self, a, b):
        return self.mult(a, self.inv(b))

    # def extended_euclidean_gcd(self, m, n):
    #     a = min(m, n)
    #     b = max(m, n)
    #     s = 0
    #     t = 1
    #     r = b
    #     old_s = 1
    #     old_t = 0
    #     old_r = a
    #     while r != 0:
    #         tmp = r
    #         quotient, r = self.euclidean_division(old_r, r)
    #         old_r = tmp
    #         old_s, s = (s, old_s - self.mult(quotient, s))
    #         old_t, t = (t, old_t - self.mult(quotient, t))
    #     return old_s, old_t, old_r
    #
    # def euclidean_division(self, m, n):
    #     """
    #     Finding q, r such that a = bq + r with deg(r) < deg(b)
    #     """
    #     a = min(m, n)
    #     b = max(m, n)
    #     q = 1
    #     r = self.add(a, self.mult(b, q))  # substraction and addition are the same in field of characteristic 2.
    #     while r.bit_length() >= b.bit_length():
    #         print(a, self.mult(b, q), r)
    #         q += 1
    #         r = self.add(a, self.mult(b, q))
    #     return q, r

    def _is_elem(self, k):
        if not (0 <= k < self.order):
            raise ValueError(f"Require 0 <= k < {self.order}, but received k = {k}.")
        return True

    def ensure_field(self):
        matrix = np.array([[self.mult(i, j) for j in range(self.order)] for i in range(self.order)])
        assert is_latin_square(matrix[1:, 1:])  # divisibility of units
        assert np.all(matrix.transpose() == matrix)  # commutativity of multiplication
        # addition if XOR, so it is commutative.
        return True


def is_latin_square(matrix):
    """
    Check that every row and columns contain the same set of elements
    and that all elements appear once.
    :param matrix:
    :return:
    """
    matrix = np.array(matrix)
    num_row, num_col = matrix.shape
    if num_row != num_col:
        return False
    n = num_col
    elems = set(matrix[0, :])
    for i in range(n):
        row_elems, row_elem_counts = np.unique(matrix[i, :], return_counts=True)
        row_elems = set(row_elems)
        if row_elems != elems or not all(row_elem_counts == 1):
            return False
    for j in range(n):
        col_elems, col_elem_counts = np.unique(matrix[:, j], return_counts=True)
        col_elems = set(col_elems)
        if col_elems != elems or not all(col_elem_counts == 1):
            return False
    return True

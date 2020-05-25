from binary_galois_field import BinaryGaloisField, is_latin_square

ord10_cayley_table = [
    [0, 3, 1, 7, 5, 9, 8, 6, 4, 2],
    [7, 0, 9, 2, 1, 5, 4, 8, 6, 3],
    [4, 2, 0, 6, 8, 7, 1, 3, 5, 9],
    [1, 7, 5, 0, 9, 8, 3, 4, 2, 6],
    [6, 1, 2, 3, 0, 4, 5, 9, 7, 8],
    [3, 6, 7, 4, 2, 0, 9, 5, 8, 1],
    [5, 8, 6, 9, 7, 2, 0, 1, 3, 4],
    [8, 9, 4, 5, 3, 6, 2, 0, 1, 7],
    [9, 4, 3, 8, 6, 1, 7, 2, 0, 5],
    [2, 5, 8, 1, 4, 3, 6, 7, 9, 0]
]


class WTAQuasigroup(object):
    def __init__(self, order=2):
        self.order = order
        if order <= 2 or order == 6:
            raise ValueError("There are no totally antisymmetric quasigroup of order 2 and 6.")
        elif order == 10:
            self.mult = self._decimal_mult
        elif order % 4 == 2:
            raise Exception("Not implemented for n = 2 mod 4")
        elif order % 2 == 1:
            self.mult = self._odd_order_mult  # -i + j mod order
        else:
            # order = 2^deg x odd_prime_prod
            self.deg = 0
            self.odd_prime_prod = order
            while self.odd_prime_prod % 2 == 0:
                self.deg += 1
                self.odd_prime_prod //= 2
            self.bgf = BinaryGaloisField(self.deg)
            if self.odd_prime_prod == 1:
                self.mult = self._binary_power_mult
            else:
                self.mult = self._composite_mult

    def inv(self, k):
        for p in range(self.order):
            if self.mult(p, k) == 0:
                return p
        return

    def _decimal_mult(self, i, j):
        return ord10_cayley_table[i][j]

    def _odd_order_mult(self, i, j):
        return (-i + j) % self.order

    def _binary_power_mult(self, i, j):
        return self.bgf.add(self.bgf.mult(2, i), j)

    def _composite_mult(self, i, j):
        # (a1, a2) * (b1, b2) = (c1, c2)
        a1, a2 = self._to_cross_product(i)
        b1, b2 = self._to_cross_product(j)
        c1 = self.bgf.add(self.bgf.mult(2, a1), b1)
        c2 = (-a2 + b2) % self.odd_prime_prod
        return c1 * self.odd_prime_prod + c2

    def _to_cross_product(self, k):
        q = k // self.odd_prime_prod
        r = k % self.odd_prime_prod
        return q, r

    def _ensure_wta_quasigroup(self):
        return


def check_weakly_totally_antisymmetric_quasigroup(matrix):
    """
    A stupid n^3 algorithm.
    But not sure if o(n^3) algorithm even exist.
    """
    if not is_latin_square:
        print("Not Latin square.")
        return False
    order = len(matrix)
    weakness = False
    for i in range(order):
        for j in range(order):
            ij = matrix[i][j]
            ji = matrix[j][i]
            if i != j and ij == ji and not weakness:
                print(f"Not strongly antisymmetric. Witness: {i}*{j} = {j}*{i} = {ij}")
                weakness = True
            for k in range(order):
                ik = matrix[i][k]
                ij_k = matrix[ij][k]
                ik_j = matrix[ik][j]
                if j != k and ij_k == ik_j:
                    print(f"Not weakly totally antisymmetric. "
                          f"Witness: ({i} * {j}) * {k} = ({i} * {k}) * {j} = {ij_k} but {j} != {k}.")
                    return False
    return True

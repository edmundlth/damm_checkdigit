import numpy as np
from wta_quasigroup import WTAQuasigroup, ord10_cayley_table, check_weakly_totally_antisymmetric_quasigroup


class Damm(object):
    def __init__(self, alphabets):
        self.alphabets = alphabets
        self._check_alphabet_set()
        self.neutral_elem = alphabets[0]  # first element is always the neutral element.
        self.order = len(alphabets)
        self.wta_quasigroup = WTAQuasigroup(self.order)
        return

    def _check_string_alphabets(self, string):
        if not set(string).issubset(set(self.alphabets)):
            raise ValueError(f"Foreign character in string detected.")

    def _check_alphabet_set(self):
        assert len(set(self.alphabets)) == len(self.alphabets), f"The set of alphabets need to be distinct."
        assert all([type(c) in [str, chr] and len(c) == 1 for c in self.alphabets]), f"All alphabets need to be of type chr."

    def _to_num_list(self, string):
        return [self._to_num(c) for c in string]

    def _to_num(self, char):
        return self.alphabets.index(char)

    def _to_char(self, num):
        return self.alphabets[num]

    def damm_process(self, string):
        self._check_string_alphabets(string)
        num_list = self._to_num_list(string)
        result = 0
        for num in num_list:
            result = self.wta_quasigroup.mult(result, num)
        return self._to_char(result)

    def damm_check(self, string):
        return self.damm_process(string) == self.neutral_elem

    def get_check_char(self, string):
        c = self.damm_process(string)
        return self._to_char(self.wta_quasigroup.inv(self._to_num(c)))

    def dammify(self, string):
        return string + self.get_check_char(string)


# def damm_process(num: str, matrix) -> int:
#     assert num.isdigit(), f"Input must be string of digit: {num}"
#     c = 0
#     for digit in num:
#         c = matrix[c][int(digit)]
#     return c
#
#
# def damm_check(num: str) -> bool:
#     return damm_process(num) == 0
#
#
# def damm_code(num: str, matrix) -> str:
#     c = damm_process(num)
#     # return inverse of c
#     possible_inverses = [a for a in matrix[c] if matrix[c][a] == 0]
#     assert len(possible_inverses) == 1, f"Inverses: {possible_inverses}"
#     check_digit = possible_inverses[0]
#     return f"{num}{check_digit}"


def large_test(num_test=1000, max_length=50, num_mutation_test=10):
    char_set = [str(a) for a in range(10)]  # digits 0 - 9
    damm = Damm(char_set)
    for _ in range(num_test):
        length = np.random.randint(1, max_length)
        test_string = ''.join([char_set[np.random.randint(0, 10)] for _ in range(length)])
        code = damm.dammify(test_string)
        org_pass_check = damm.damm_check(code)
        assert org_pass_check, f"Org string={test_string}\tCode={code}"
        for _ in range(num_mutation_test):
            mutated = list(code)
            rand_pos = np.random.randint(0, length)
            rand_char = char_set[np.random.randint(0, 10)]
            while rand_char == mutated[rand_pos]:
                rand_char = char_set[np.random.randint(0, 10)]
            mutated[rand_pos] = rand_char
            mutated = ''.join(mutated)
            assert not damm.damm_check(mutated), f"Org code={code}\nMutated={mutated}"

            swapped = list(code)
            if rand_pos == length - 1:
                rand_pos -= 1
            if swapped[rand_pos] != swapped[rand_pos + 1]:
                tmp = swapped[rand_pos]
                swapped[rand_pos] = swapped[rand_pos + 1]
                swapped[rand_pos + 1] = tmp
                swapped = ''.join(swapped)
                assert not damm.damm_check(swapped), f"Org code={code}\nSwapped={swapped}"
    return


if __name__ == "__main__":
    test_string = "0123456789"
    print([type(a) for a in list(set(test_string))])

    damm = Damm(list(set(test_string)))
    code = damm.dammify(test_string)
    org_pass_check = damm.damm_check(code)
    mutated_codes = [
        "1023456789" + code[-1],
        "0023456789" + code[-1]
    ]
    print(f"Original digit string is {test_string}")
    print(f"which has Damm code {code}")
    print(f"It pass check digit test? {org_pass_check}\n")
    for mutated in mutated_codes:
        pass_test = damm.damm_check(mutated)
        print(f"Now mutate to {mutated}")
        print(f"does it still pass the test? {pass_test}")
        assert not pass_test
    assert org_pass_check
    assert check_weakly_totally_antisymmetric_quasigroup(ord10_cayley_table)
    large_test()

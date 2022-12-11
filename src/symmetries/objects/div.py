from copy import deepcopy


class Div():

    def __init__(self, numerator, denominator) -> None:
        self.numerator = numerator
        self.denominator = denominator

    def __eq__(self, other):
        if isinstance(other, Div):
            num = self.numerator == other.numerator
            den = self.denominator == other.denominator

            return (num and den)

        else:
            return False

    def __repr__(self):
        num = ''
        try:
            assert self.denominator.terms
            if len(self.denominator.terms) > 1 or self.denominator.coefficient != 1:
                num = f'({self.denominator})'
        except AttributeError:
            pass
        if not num:
            num = f'{self.denominator}'
        return f'{self.numerator}/{num}'

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            try:
                assert self.denominator.coefficient
                self_copy = deepcopy(self)
                self_copy.denominator.coefficient = 1
                self_copy.numerator *= other / self.denominator.coefficient
                return self_copy

            except AttributeError:
                self_copy = deepcopy(self)
                self_copy.numerator *= other
                return self_copy
        else:
            return other*self

    def __rmul__(self, other):
        return self*other

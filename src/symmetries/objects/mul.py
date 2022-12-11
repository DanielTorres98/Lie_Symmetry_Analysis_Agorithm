from copy import deepcopy
from .add import Add
from .div import Div


class Mul():
    """Class for multiplication of multiple terms."""

    def __init__(self, terms: tuple, coefficient: float = 1) -> None:
        self.coefficient = coefficient
        # Tuple[Union[Add, Mul, Variable, DependentVariable]]
        self.terms = terms

    def __eq__(self, other):
        if isinstance(other, Mul):
            return self.terms == other.terms
            # they have to be in the same order, need a way to sort them
        else:
            if len(self.terms)==1:
                return other in self.terms
            else:
                return False

    def __repr__(self):
        display = f'{round(self.coefficient,4)}' if self.coefficient != 1 else ''
        if display == '-1':
            display = '-'
        for term in self.terms:
            display += f'{term.__repr__()}*'
        return display[:-1]

    def __mul__(self, other):
        """Multiplication method for variables."""
        if isinstance(other, (int, float)):
            if other:
                res = deepcopy(self)
                res.coefficient *= other
                return res
            else:
                return 0

        elif isinstance(other, Mul):
            coefficient = self.coefficient*other.coefficient
            terms = self.terms+other.terms
            base = terms[0]
            for term in terms[1:]:
                base = base*term
            if isinstance(base, Mul):
                return Mul(base.terms, coefficient)
            else:
                return Mul((base,), coefficient)

        elif isinstance(other, Add):
            addition_terms = []
            for term in other.terms:
                addition_terms.append(self*term)
            return Add(tuple(addition_terms))

        elif isinstance(other, Div):
            if isinstance(other.denominator, Mul):
                if self in other.denominator.terms:
                    other_copy = deepcopy(other)
                    other.denominator.terms.remove(self)

                    other_copy.denominator.terms = tuple(other.denominator.terms)
                    other_copy.numerator *= other.denominator.coefficient
                    return other_copy
            else:
                other_copy = deepcopy(other)
                other_copy.numerator = other_copy.numerator * self
                return other_copy


        else: # variable
            result = deepcopy(self)
            if other in result.terms:
                terms = []
                for term in result.terms:
                    if term==other:
                        terms.append(term*other)
                    else:
                        terms.append(term)
                return Mul(tuple(terms), self.coefficient)
            else:
                terms = list(result.terms)+[other]
                terms.sort(key=lambda x: x.name)
                result.terms = tuple(terms)
                return result

    def __rmul__(self, other):
        return self*other

    def __radd__(self, other):
        return self+other

    def __add__(self, other):
        if isinstance(other, (int, float)):
            if other:
                return Add((self, other))
            else:
                return self

        elif isinstance(other, Mul):
            if other == self:
                res = deepcopy(self)
                res.coefficient = self.coefficient + other.coefficient
                if res.coefficient:
                    return res
                else:
                    return 0
            else:
                return Add((self, other))

        elif isinstance(other, Add):
            addition_terms = []
            accounted_for_self = False
            for term in other.terms:
                if isinstance(term, Mul) and self == term:
                    addition_terms.append(term+self)
                    accounted_for_self = True
                elif term == self:
                    addition_terms.append(term+self)
                    accounted_for_self = True
                else:
                    addition_terms.append(term)
            if not accounted_for_self:
                addition_terms.append(self)
            return Add(tuple(addition_terms))
        
        else:
            if other == self:
                coeff = self.coefficient+1
                return Mul(self.terms, coeff)
            else:
                return Add((self, other))

    def __rsub__(self, other):
        result = deepcopy(self)
        result.coefficient *= -1
        return other+result

    def __sub__(self, other):
        return self+(-1*other)

    def __pow__(self, other):
        assert isinstance(other, int), 'No partial power.'

        if other==0:
            return 1

        elif other>0:
            result = deepcopy(self)
            result.coefficient = result.coefficient ** other
            terms = []
            for term in result.terms:
                terms.append(term**other)

            result.terms = tuple(terms)
            return result

        else: # power<0
            power_mag = abs(other)
            result = self**power_mag
            return Div(1, result)

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            return Div(other, self)  # TODO: change this
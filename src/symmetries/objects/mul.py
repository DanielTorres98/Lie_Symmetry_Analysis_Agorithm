from variables import Variable
from copy import deepcopy
from add import Add


class Mul():
    """Class for multiplication of multiple terms."""

    def __init__(self, terms: tuple = (), coefficient: float = 1) -> None:
        self.coefficient = coefficient
        # Tuple[Union[Add, Mul, Variable, DependentVariable]]
        self.terms = terms

    def __eq__(self, other):
        if isinstance(other, Mul):
            return (self.terms == other.terms)
            # they have to be in the same order, need a way to sort them
        elif isinstance(other, Variable):
            return other == self
        else:
            return False

    def __repr__(self):
        display = f'{self.coefficient}' if self.coefficient != 1 else ''
        if display == '-1':
            display = '-'
        for term in self.terms:
            display += f'{term.__repr__()}*'
        return display[:-1]

    def __mul__(self, other):
        """Multiplication method for variables."""
        if isinstance(other, Variable):
            terms = []
            inside = False
            for term in self.terms:
                if term.name == other.name:
                    terms.append(other*term)
                    inside = True
                else:
                    terms.append(term)
            if not inside:
                terms.append(other)
            res = deepcopy(self)
            terms.sort(key=lambda x: x.name)
            res.terms = tuple(terms)
            return res

        elif isinstance(other, float) or isinstance(other, int):
            if other:
                res = deepcopy(self)
                res.coefficient *= other
                return res
            else:
                return 0

        elif isinstance(other, Mul):
            coeff = self.coefficient*other.coefficient
            terms = self.terms+other.terms
            base = terms[0]
            for term in terms[1:]:
                base = base*term
            if coeff:
                return Mul(base.terms, coefficient=coeff)
            else:
                0

        elif isinstance(other, Add):
            addition_terms = []
            for term in other.terms:
                addition_terms.append(self*term)
            return Add(tuple(addition_terms))

    def __rmul__(self, other):
        return self*other

    def __add__(self, other):
        if isinstance(other, Variable):
            if other == self:
                coeff = self.coefficient+1
                return Mul((self.terms[0],), coeff)
            else:
                return Add((self, other))

        elif isinstance(other, float) or isinstance(other, int):
            if other:
                return Add((self, other))
            else:
                return self

        elif isinstance(other, Mul):
            if other == self:
                res = deepcopy(self)
                res.coefficient += other.coefficient
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
                elif isinstance(term, Variable):
                    if term == self:
                        addition_terms.append(term+self)
                        accounted_for_self = True
                else:
                    addition_terms.append(term)
            if not accounted_for_self:
                addition_terms.append(self)
            return Add(tuple(addition_terms))

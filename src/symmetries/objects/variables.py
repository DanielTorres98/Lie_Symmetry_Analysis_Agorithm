from symmetries.utils.constants import greek_alphabet
from .add import Add
from .mul import Mul
from copy import deepcopy


class Variable():
    """Class for independent variables."""
    greek_alphabet = greek_alphabet

    def __init__(self, name: str) -> None:
        self.name = name

        self.symbol = self.name
        if self.name in self.greek_alphabet:
            self.symbol = self.greek_alphabet[self.name]

    def __repr__(self):
        return self.symbol

    def __eq__(self, other):
        if isinstance(other, Variable):
            return (self.name == other.name)
        elif isinstance(other, Mul):
            if len(other.terms) == 1:
                return other.terms[0] == self
        else:
            return False

    def __mul__(self, other):
        """Multiplication method for variables."""
        if isinstance(other, Variable):
            if other.name == self.name:
                if isinstance(other, Power):
                    power = other.power
                else:
                    power = 1
                if isinstance(self, Power):
                    power += self.power
                    return Power(self.term, power)
                else:
                    power += 1
                return Power(self, power)

            else:
                terms = [self, other]
                terms.sort(key=lambda x: x.name)
                return Mul(tuple(terms))

        if isinstance(other, (int, float)):
            if other:
                return Mul((self,), coefficient=other)
            else:
                return 0

        elif isinstance(other, Mul):
            res = deepcopy(self)
            for term in other.terms:
                res = res*term

            if isinstance(res, Variable):
                return Mul(terms=(res,), coefficient=other.coefficient)
            else:
                terms = list(res.terms)
                terms.sort(key=lambda x: x.name)
                return Mul(terms=tuple(terms), coefficient=other.coefficient)

        elif isinstance(other, Add):
            results = []
            for term in other.terms:
                res = self*term
                results.append(res)
            return Add(terms=tuple(results))

    def __rmul__(self, other):
        return self*other

    def __pow__(self, other: int):
        if isinstance(self, Power):
            power = self.power*other
            return Power(self.term, power)
        else:
            return Power(self, other)

    def __sub__(self, other):
        return self+(-1*other)

    def __rsub__(self, other):
        return -1*self+other

    def __radd__(self, other):
        return self+other

    def __add__(self, other):
        """Addition method for variables."""
        if isinstance(other, Variable):
            if other == self:
                return Mul((self,), coefficient=2)
            else:
                return Add((self, other))

        elif isinstance(other, (int, float)):
            if other:
                return Add((self, other))
            else:
                return self

        elif isinstance(other, Mul):
            if len(other.terms) == 1 and other.terms[0] == self:
                coefficient = other.coefficient + 1
                if coefficient:
                    return Mul((self,), coefficient)
                else:
                    return 0
            else:
                return Add((self, other))

        elif isinstance(other, Add):
            return other+self


class Power(Variable):

    def __init__(self, term, power):
        self.power = power
        self.term = term
        self.name = term.name

    def __eq__(self, other):
        if isinstance(other, Mul):
            if len(other.terms) == 1:
                return other.terms[0] == self
        elif isinstance(other, Power):
            return (self.name == other.name and self.power == other.power)
        else:
            return False

    def __repr__(self):
        display = self.term.__repr__()
        if self.power and self.power != 1:
            display += f'^{self.power}'
        return display


class Function(Variable):
    """Class for functions, child class of variables."""

    def __init__(self, name: str, dependencies: tuple = (), derivatives: tuple = ()) -> None:
        self.derivatives = derivatives
        self.dependencies = dependencies

        super().__init__(name)

    def __eq__(self, other):
        if isinstance(other, Function):
            return (self.name == other.name and self.derivatives == other.derivatives)
        else:
            return False

    def __repr__(self):
        display = f'{self.symbol}('
        for term in self.dependencies:
            display += f'{term},'
        display = display[:-1]+')'
        for term in self.derivatives:
            display += f'_{term}'
        return display

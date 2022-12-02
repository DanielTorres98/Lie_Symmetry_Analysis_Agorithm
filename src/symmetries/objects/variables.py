from symmetries.utils.constants import greek_alphabet
from typing import List, Union
from copy import deepcopy

class Variable():
    """Class for independent variables."""
    greek_alphabet = greek_alphabet

    def __init__(self, name: str, power: int = 1) -> None:
        self.name = name
        
        self.symbol = self.name
        if self.name in self.greek_alphabet:
            self.symbol = self.greek_alphabet[self.name]
        self.power = power

    def __repr__(self):
        display = self.symbol
        if self.power and self.power != 1:
            display += f'^{self.power}'
        return display

    def __eq__(self, other):
        return (self.name == other.name)

    def __mul__(self, other):
        """Multiplication method for variables."""
        if isinstance(other, DependentVariable):
            if other.name == self.name:
                res = deepcopy(self)
                res.power += other.power
                return res
            else:
                return Mul((self, other))

        elif isinstance(other, Variable):
            if other.name == self.name:
                return Variable(self.name, other.power+self.power)
            return Mul((self, other))

        elif isinstance(other, float) or isinstance(other, int):
            return Mul((self,), coeff=other)

        elif isinstance(other, Mul):
            results = []
            for term in other.terms:
                res = self.__mul__(term)
                results.append(res)
            return Mul(terms=tuple(results), coefficient=other.coefficient)

        elif isinstance(other, Add):
            results = []
            for term in other.terms:
                res = self.__mul__(term)
                results.append(res)

            return Add(terms=tuple(results))


class DependentVariable(Variable):
    """Class for dependent variables, child class of variables."""

    def __init__(self, name: str, power: int = 1, dependencies: tuple = (), derivatives: tuple = ()) -> None:
        self.derivatives = derivatives
        self.dependencies = dependencies
        
        self.symbol = name
        if name in self.greek_alphabet:
            self.symbol = self.greek_alphabet[self.name]
        
        self.power = power
        self.name = self.__repr__().split('^')[0]


    def __eq__(self, other):
        return (self.name == other.name and self.derivatives == other.derivatives)

    def __repr__(self):
        display = f'{self.symbol}('
        for term in self.dependencies:
            display += f'{term},'
        display = display[:-1]+')'
        if self.power and self.power != 1:
            display += f'^{self.power}'
        for term in self.derivatives:
            display += f'_{term}'
        return display


class Mul():
    """Class for multiplication of multiple terms."""

    def __init__(self, terms: tuple = (), coefficient: float = 1) -> None:
        self.coefficient = coefficient
        # Tuple[Union[Add, Mul, Variable, DependentVariable]]
        self.terms = terms

    def __repr__(self):
        display = f'{self.coefficient}' if self.coefficient != 1 else ''
        for term in self.terms:
            display += f'{term.__repr__()}*'
        return display[:-1]

    def __mul__(self, other):
        """Multiplication method for variables."""

        if isinstance(other, Variable):
            terms = []
            for term in self.terms:
                if term.name==other.name:
                    terms.append(other*term)
                else:
                    terms.append(term)
            res = deepcopy(self)
            res.terms = tuple(terms)
            return res

        elif isinstance(other, float) or isinstance(other, int):
            res = deepcopy(self)
            res.coefficient *= other
            return res

        elif isinstance(other, Mul):
            pass

        elif isinstance(other, Add):
            pass
class Add():
    """Class for addition of multiple terms."""

    def __init__(self, terms: tuple = ()) -> None:
        self.terms = terms

    def __repr__(self):
        display = ''
        for term in self.terms:
            display += f'{term.__repr__()}+'
        return display[:-1]

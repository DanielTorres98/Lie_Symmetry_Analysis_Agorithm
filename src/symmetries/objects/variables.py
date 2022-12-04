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
                res = deepcopy(self)
                res.power += other.power
                return res
            else:
                return Mul((self, other))

        elif isinstance(other, float) or isinstance(other, int):
            return Mul((self,), coefficient=other)

        elif isinstance(other, Mul):
            results = []
            for term in other.terms:
                res = self*term
                results.append(res)
            return Mul(terms=tuple(results), coefficient=other.coefficient)

        elif isinstance(other, Add):
            results = []
            for term in other.terms:
                res = self*term
                results.append(res)
            return Add(terms=tuple(results))

    def __add__(self, other):
        """Addition method for variables."""
        if isinstance(other, Variable):
            if other == self:
                return Mul((self,), coefficient=2)
            else:
                return Add((self, other))

        elif isinstance(other, float) or isinstance(other, int):
            return Add((self, other))

        elif isinstance(other, Mul):
            if len(other.terms) == 1 and other.terms[0] == self:
                coeff = other.coefficient + 1
                return Mul((self,), coefficient=coeff)
            return Add((self, other))

        elif isinstance(other, Add):
            return Add(terms=other.terms+(self,))


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
        if isinstance(other, DependentVariable):
            return (self.name == other.name and self.derivatives == other.derivatives and self.power == other.power)
        else:
            return False

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

    def __eq__(self, other):
        if isinstance(other, Mul):
            return (self.terms == other.terms)
            # they have to be in the same order, need a way to sort them
        elif isinstance(other, Variable):
            return other == self
        else:
            return False

    def __repr__(self):
        if self.coefficient == 0:
            return '0'

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
            res.terms = tuple(terms)
            return res

        elif isinstance(other, float) or isinstance(other, int):
            res = deepcopy(self)
            res.coefficient *= other
            return res

        elif isinstance(other, Mul):
            coeff = self.coefficient*other.coefficient
            terms = []
            for term_outer in self.terms:
                for term_inner in other.terms:
                    terms.append(term_outer*term_inner)
            return Mul(tuple(terms), coefficient=coeff)

        elif isinstance(other, Add):
            addition_terms = []
            for term in other.terms:
                addition_terms.append(self*term)
            return Add(tuple(addition_terms))

    def __add__(self, other):
        if isinstance(other, Variable):
            if other == self:
                coeff = self.coefficient+1
                return Mul((self.terms[0],), coeff)
            else:
                return Add(terms=(self, other))

        elif isinstance(other, float) or isinstance(other, int):
            return Add(terms=(self, other))

        elif isinstance(other, Mul):
            if other == self:
                res = deepcopy(self)
                res.coefficient += other.coefficient
                return res
            else:
                return Add(terms=(self, other))

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


class Add():
    """Class for addition of multiple terms."""

    def __init__(self, terms: tuple = ()) -> None:
        self.terms = terms

    def __repr__(self):
        display = ''
        for term in self.terms:
            display += f'{term.__repr__()}+'
        return display[:-1]

    def __add__(self, other):
        if isinstance(other, Add):
            res = deepcopy(self)
            for term in other.terms:
                res = res+term
            return res

        else:
            if other in self.terms:
                addition_terms = [term if term !=
                                  other else term+other for term in self.terms]
                return Add(tuple(addition_terms))
            else:
                return Add(self.terms+(other,))

    def __mul__(self, other):
        if isinstance(other, Add):
            additions = [self*term for term in other.terms]
            base = additions[0]
            for expansion in additions[1:]:
                base = base+expansion
            return base

        elif isinstance(other, float) or isinstance(other, int):
            return Add(tuple([term*other for term in self.terms]))
        else:
            return Add(tuple([other*term for term in self.terms]))

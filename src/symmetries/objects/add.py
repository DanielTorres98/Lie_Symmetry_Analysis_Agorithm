from copy import deepcopy

class Add():
    """Class for addition of multiple terms."""

    def __init__(self, terms: tuple = ()) -> None:
        self.terms = terms

    def __repr__(self):
        display = ''
        for term in self.terms:
            if term.__repr__()[0]=='-':
                pass
            else:
                if display:
                    display += '+'
            display += f'{term.__repr__()}'
        return display

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

    def __rsub__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return self-other

    def __sub__(self, other):
        if isinstance(other, Add):
            res = deepcopy(self)
            for term in other.terms:
                res = res-term
            return res
        else:
            if other in self.terms:
                addition_terms = [term if term !=
                                  other else term-other for term in self.terms]
                return Add(tuple(addition_terms))
            else:
                term = 0-other
                return Add(self.terms+(term,))

    def __mul__(self, other):
        if isinstance(other, Add):
            additions = []
            for inner_term in self.terms:
                for outer_term in other.terms:
                    additions.append(inner_term*outer_term)

            base = additions[0]
            for expansion in additions[1:]:
                base = base+expansion
            return base

        elif isinstance(other, float) or isinstance(other, int):
            if other:
                return Add(tuple([term*other for term in self.terms]))
            else:
                return 0
        else:
            terms = [term*other for term in self.terms]
            return Add(tuple(terms))

    def __rmul__(self, other):
        return self*other
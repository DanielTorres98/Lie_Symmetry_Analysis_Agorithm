from copy import deepcopy


class Add():
    """Class for addition of multiple terms."""

    def __init__(self, terms: tuple = ()) -> None:
        self.terms = terms

    def __repr__(self):
        display = ''
        for term in self.terms:
            if term.__repr__()[0] == '-':
                pass
            else:
                if display:
                    display += '+'
            display += f'{term.__repr__()}'
        return display

    def __radd__(self, other):
        return self+other

    def __add__(self, other):
        if isinstance(other, Add):
            res = deepcopy(self)
            for term in other.terms:
                res = res+term
            return res

        else:
            if isinstance(other, (int, float)):
                if any(isinstance(term, (int, float)) for term in self.terms):
                    addition_terms = [
                        term + other if isinstance(term, (int, float)
                                                ) else term for term in self.terms]
                    return Add(tuple(addition_terms))

            elif other in self.terms:
                addition_terms  = []
                for term in self.terms:
                    if term != other:
                        addition_terms.append(term)
                    else:
                        res = term+other
                        if res:
                            addition_terms.append(res)

                return Add(tuple(addition_terms))
            
            return Add(self.terms+(other,))

    def __rsub__(self, other):
        return -1*self+other

    def __sub__(self, other):
        return self+(-1*other)

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

        elif isinstance(other, (int, float)):
            if other:
                return Add(tuple([term*other for term in self.terms]))
            else:
                return 0
        else:
            terms = [term*other for term in self.terms]
            return Add(tuple(terms))

    def __rmul__(self, other):
        return self*other

    def __pow__(self, other):
        assert isinstance(other, int)
        self_copy = deepcopy(self)

        for _ in range(other-1):
            self_copy = self_copy * self

        return self_copy

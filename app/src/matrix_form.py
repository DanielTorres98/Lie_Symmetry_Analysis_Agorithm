
"""Module containing functions for turning the equations into matrix form"""
import numpy as np
import sympy as sp
from copy import deepcopy

class MatrixForm():

    def __init__(self, latex_form) -> None:
        self.equations = self._split(latex_form)
        self.variables = self._get_variables()
        self.matrix = sp.zeros(len(self.equations), len(self.variables))

    def _split(self, latex_form):
        equations = latex_form.split("\\\n")
        equations = [eq.split('=')[0] for eq in equations]
        return [eq.split('+') for eq in equations if eq]

    def _get_variables(self):
        variables = []
        for equation_terms in self.equations:
            for term in equation_terms:
                variable = term.split('*')[-1]
                variable_s = sp.symbols(variable)
                if variable_s not in variables:
                    variables.append(variable_s)
        return variables
    
    def populate_matrix(self):
        self.matrix = sp.zeros(len(self.equations), len(self.variables))
        for r, equation_terms in enumerate(self.equations):
            for term in equation_terms:
                coefficients = term.split('*')[:-1]
                if len(equation_terms)>1:
                    coefficient = int(coefficients[0])
                    if len(coefficients)>1:
                        coefficient *= sp.symbols(coefficients[1])
                else:
                    coefficient = 1
                idx = self.variables.index(sp.symbols(term.split('*')[-1]))
                if not self.matrix[r,idx]:
                    self.matrix[r,idx] = coefficient
                else:
                    self.matrix[r,idx] += coefficient

    def identify_common_terms(self):
        terms = []
        non_zero = []
        for i in range(self.matrix.shape[0]):
            terms.append(np.count_nonzero(self.matrix[i,:]))
            non_zero.append(list(np.nonzero(self.matrix[i,:])[1]))
        zipped = zip(terms, range(self.matrix.shape[0]), non_zero) 
        return sorted(zipped, key = lambda x: x[0])
    
    def delete_single_terms(self, inplace=True):
        m = deepcopy(self.matrix)
        v = deepcopy(self.variables)

        i = 0
        for n, r, c in self.identify_common_terms():
            if n>1:
                break
            m.row_del(r-i)
            m.col_del(c[0]-i)
            del v[c[0]-i]
            i += 1
        if inplace:
            self.matrix = m
            self.variables = v
        else:
            return m, v
        
# latex_form = latex.format_determining_equations(False)
# m = MatrixForm(latex_det)
# m.populate_matrix()
# m.identify_common_terms()
# m.matrix[13,:] = m.matrix[13,:]+m.matrix[15,:]-m.matrix[19,:]*(sp.symbols('\\rho')-sp.symbols('b'))
# m.matrix.row_del(19)
# m.matrix.row_del(15)
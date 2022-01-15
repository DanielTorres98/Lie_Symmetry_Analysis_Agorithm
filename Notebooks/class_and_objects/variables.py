import numpy as np
from copy import deepcopy


class Variable():
    greek_alphabet = {
    'Alpha'  : u'\u0391',
    'Beta'   : u'\u0392',
    'Gamma'  : u'\u0393',
    'Delta'  : u'\u0394',
    'Epsilon': u'\u0395',
    'Zeta'   : u'\u0396',
    'Eta'    : u'\u0397',
    'Theta'  : u'\u0398',
    'Iota'   : u'\u0399',
    'Kappa'  : u'\u039A',
    'Lamda'  : u'\u039B',
    'Mu'     : u'\u039C',
    'Nu'     : u'\u039D',
    'Xi'     : u'\u039E',
    'Omicron': u'\u039F',
    'Pi'     : u'\u03A0',
    'Rho'    : u'\u03A1',
    'Sigma'  : u'\u03A3',
    'Tau'    : u'\u03A4',
    'Upsilon': u'\u03A5',
    'Phi'    : u'\u03A6',
    'Chi'    : u'\u03A7',
    'Psi'    : u'\u03A8',
    'Omega'  : u'\u03A9',
    'alpha'  : u'\u03B1',
    'beta'   : u'\u03B2',
    'gamma'  : u'\u03B3',
    'delta'  : u'\u03B4',
    'epsilon': u'\u03B5',
    'zeta'   : u'\u03B6',
    'eta'    : u'\u03B7',
    'theta'  : u'\u03B8',
    'iota'   : u'\u03B9',
    'kappa'  : u'\u03BA',
    'lamda'  : u'\u03BB',
    'mu'     : u'\u03BC',
    'nu'     : u'\u03BD',
    'xi'     : u'\u03BE',
    'omicron': u'\u03BF',
    'pi'     : u'\u03C0',
    'rho'    : u'\u03C1',
    'sigma'  : u'\u03C3',
    'tau'    : u'\u03C4',
    'upsilon': u'\u03C5',
    'phi'    : u'\u03C6',
    'chi'    : u'\u03C7',
    'psi'    : u'\u03C8',
    'omega'  : u'\u03C9',
}

    def __init__(self, name, func=False, derivatives={}):
        self.name = name
        self.symbol = self.name
        if self.name in self.__class__.greek_alphabet:
            self.symbol = self.__class__.greek_alphabet[self.name]
        self.func = func
        self.derivatives = derivatives

    def __str__(self):
        display = self.symbol
        if not self.derivatives:
            return display
        else:
            for term, times in self.derivatives.items():
                for n in range(times):
                    display += f'_{term}'
            return display
    
    def __repr__(self):
        display = self.symbol
        if not self.derivatives:
            return display
        else:
            for term, times in self.derivatives.items():
                for n in range(times):
                    display += f'_{term}'
            return display
    
    def __eq__(self, other):
        return (self.name == other.name and self.derivatives ==  other.derivatives)
    
    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return term(float(other), [self])
        elif isinstance(other, Variable):
            if self == other:
                return term(1.0, [self], exponents=[2])
            else:
                return term(1.0, [self, other])

    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return term(float(other), [self])
        elif isinstance(other, Variable):
            if self == other:
                return term(1.0, [self], exponents=[2])
    
    def __add__(self, other):
        pass
            
    
    def D(self, *others):
        # Derivative method
        for other in others:
            if not isinstance(other, Variable):
                raise Exception(f'Not possible to perform with respect to {other} \n'
                                f'Please provide an object Variable as input.') 
        if self.func:
            derivatives = deepcopy(self.derivatives)
            if self in others:
                if len(others) > 1:
                    return 0
                else: 
                    return 1
            for other in others: 
                if other.symbol not in derivatives:
                    derivatives[other.symbol] = 1
                else:
                    derivatives[other.symbol] += 1 
            return Variable(self.symbol, True, derivatives)
        else:
            if self in others and len(others)==1:
                return 1
            return 0
            

class term():
    def __init__(self, coefficient, variables, exponents=None):
        self.coefficient = coefficient
        self.variables = variables
        if not exponents:
            self.exponents = np.ones(len(variables))
        else:
            self.exponents = np.array(exponents)
        if len(self.variables) != len(self.exponents):
            raise Exception('There should be as many exponents as variables')

    def __str__(self):
        display = ''
        if self.coefficient != 1:
            display = f'{self.coefficient}'
        vars_exps = zip(self.variables, self.exponents)
        for var, exp in vars_exps:
            display += f'{var}'
            if exp != 1:
                display += f'^{exp}'
        return display

    def __repr__(self):
        display = ''
        if self.coefficient != 1:
            display = f'{self.coefficient}'
        vars_exps = zip(self.variables, self.exponents)
        for var, exp in vars_exps:
            display += f'{var}'
            if exp != 1:
                display += f'^{exp}'
        return display

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return term(self.coefficient*other, self.variables, self.exponents)
        elif isinstance(other, Variable):
            variables = np.array(deepcopy(self.variables))
            exponents = deepcopy(self.exponents)
            if other not in variables:
                variables = variables.append(other)
                exponents = np.concatenate((self.exponents, 1))
            else:
                exponents[np.where(np.array(variables)==other)] += 1
            return term(self.coefficient, variables, exponents=exponents)
        elif isinstance(other, term):
            coefficient =  self.coefficient * other.coefficient
            variables = deepcopy(self.variables)
            exponents = deepcopy(self.exponents)
            other_variables = deepcopy(other.variables)
            other_exponents = deepcopy(other.exponents)
            idx_list = []
            for var in variables:
                idx = np.where(other.variables==var)
                if idx:
                    exponents[np.where(variables==var)] += other.exponents[idx]
                    idx_list.append(idx)
            all_vars = variables + other_variables[np.where(other_variables!=variables)]
            all_expo = exponents + other_exponents[np.where(other_exponents!=exponents)]
            return term(coefficient, all_vars, all_expo)

    def __rmul__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return term(self.coefficient*other, self.variables, self.exponents)
        elif isinstance(other, Variable):
            variables = np.array(deepcopy(self.variables))
            exponents = deepcopy(self.exponents)
            if other not in variables:
                variables = variables.append(other)
                exponents = np.concatenate((self.exponents, 1))
            else:
                exponents[np.where(np.array(variables)==other)] += 1
            return term(self.coefficient, variables, exponents=exponents)
        elif isinstance(other, term):
            coefficient =  self.coefficient * other.coefficient
            variables = deepcopy(self.variables)
            exponents = deepcopy(self.exponents)
            other_variables = deepcopy(other.variables)
            other_exponents = deepcopy(other.exponents)
            idx_list = []
            for var in variables:
                idx = np.where(other.variables==var)
                if idx:
                    exponents[np.where(variables==var)] += other.exponents[idx]
                    idx_list.append(idx)
            all_vars = variables + other_variables[np.where(other_variables!=variables)]
            all_expo = exponents + other_exponents[np.where(other_exponents!=exponents)]
            return term(coefficient, all_vars, all_expo)
            
                    
x = Variable('x')
y = Variable('y')

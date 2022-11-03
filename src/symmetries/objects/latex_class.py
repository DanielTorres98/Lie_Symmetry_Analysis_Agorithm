class Latex():
    def __init__(self, determining_equations:dict, independent_variables:list,
                 dependent_variables:list, constants:list):
        self.equations = determining_equations
        self.variables = independent_variables
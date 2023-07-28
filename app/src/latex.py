"""Module containing functions for printing equations in Latex form"""


class Latex():
    """Latex class"""

    def __init__(self, determining_equations, var_dict: dict, variables: list,
                 constants: list, replace_plus_minus: bool = True, print_ones: bool = False):
        """Constructor

        Parameters
        ----------
        determining_equations : dict
            _description_
        var_dict : dict
            dictionary translating the variable name to latex format
        variables : list
            list with all variables (independent and dependant) written in latex format
        constants : list
            list containing all constants in string format
        """
        self.determining_equations = determining_equations
        self.var_dict = var_dict
        self.variables = variables
        self.constants = constants

        self.replace_plus_minus = replace_plus_minus
        self.print_ones = print_ones

        diff = len(self.constants) - len(self.variables)
        alpha_id = [f'alpha_{idx}' for idx in range(diff)]
        self.symbolic_constants = constants + self.variables + alpha_id

    def format_equation_term(self, term: dict, one_term: bool):
        """Given a dictionary it returns the symbolic equivalent. It drops all constants
        if it is just one term.

        Parameters
        ----------
        term : dict
            dictionary containing the information of the term
        one_term : bool
            True if it is the only term in the equation

        Returns
        -------
        str
            latex code to write the term
        """
        variable = self.var_dict[term['variable']]
        derivatives = term['derivatives']

        coefficient = str(term['coefficient'])
        if self.print_ones:
            coefficient += '*'
        else:
            if coefficient == '-1':
                coefficient = "-"
            elif coefficient == '1':
                coefficient = ""

        term_latex = ''
        if not one_term:
            for cte, n in zip(self.symbolic_constants, term['constants']):
                cte = str(cte)
                if len(cte) > 1:
                    cte = "\\" + cte

                if n == 1:
                    term_latex += cte
                elif n > 1:
                    term_latex += cte + '^' + str(n)

            if term_latex:
                term_latex += '*'

        D = self.format_derivatives(derivatives, variable)
        return coefficient  + term_latex + D

    def format_derivatives(self, derivatives: list, variable: str):
        """Given a list of derivatives executes all the derivatives on the variable.

        Parameters
        ----------
        derivatives : list
            list of integers containing the order of the derivative with respect to the
            variable
        variable : str
            variable to be differentiated

        Returns
        -------
        str
            latex code for the derivatives
        """
        var_str = '\\' + variable
        for D, var in zip(derivatives, self.variables):
            var = str(var)
            for _ in range(D):
                if len(var) > 1:
                    var = "\\" + var

                if '_' not in var_str:
                    var_str += '_'
                var_str += '{' + var + '}'
        return var_str

    def format_equation(self, equation: dict) -> str:
        """This code generates the latex string format of the determinant equations

        Parameters
        ----------
        equation : dict
            dictionary containing a determinant equation

        Returns
        -------
        str
            string with all terms of a single equation typed in latex format
        """
        latex_equation = ''
        one_term = (len(equation) == 1)
        for term in equation:
            latex_equation += self.format_equation_term(term, one_term) + '+'
        return latex_equation[:-1] + '=0'

    def format_determining_equations(self) -> str:
        """Gives the latex code version of remaining determining equation.
        Returns
        -------
        str
            single string with all the equations typed in latex format
        """
        latex_system_of_eqns = ''
        for equation in self.determining_equations.values():
            latex_system_of_eqns += self.format_equation(
                equation) + "\\" + "\\" + "\n"
        if self.replace_plus_minus:
            latex_system_of_eqns = latex_system_of_eqns.replace("+-", "-")
        return latex_system_of_eqns

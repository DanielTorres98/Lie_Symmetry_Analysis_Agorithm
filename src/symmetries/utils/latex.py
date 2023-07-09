"""Module containing functions for printing equations in Latex form"""


class Latex():

    def __init__(self, determining_equations, var_dict: dict, variables: list, constants: list):
        """_summary_

        Parameters
        ----------
        determining_equations : _type_
            _description_
        var_dict : dict
            dictionary translating the variable name to latex format
        variables : list
            list with all variables (independent and dependant) with the constants
            written in latex format
        constants : list
            list containing all constants
        """
        self.determining_equations = determining_equations
        self.var_dict = var_dict
        self.variables = variables
        self.constants = constants

    def dict_to_latex(self, term: dict, symbolic_constants: list, one_term: bool):
        """Given a dictionary it returns the symbolic equivalent. It drops all constants
        if it is just one term.

        Parameters
        ----------
        term : dict
            dictionary containing the information of the term
        symbolic_constants : list
            [description]
        one_term : bool
            True if it is the only term in the equation

        Returns
        -------
        str
            The latex code to write the term. 
        """
        variable = self.var_dict[term['variable']]
        derivatives = term['derivatives']
        a = ''
        if one_term:
            coeff = ''
        else:
            for cte, n in zip(symbolic_constants, term['constants']):
                if len(str(cte)) > 1:
                    if n == 1:
                        a += "\\" + str(cte)
                    if n > 1:
                        a += "\\" + str(cte) + '^' + str(n)
                else:
                    if n == 1:
                        a += str(cte)
                    if n > 1:
                        a += str(cte) + '^' + str(n)
            coeff = str(term['coefficient'])
        D = self.format_derivatives(derivatives, variable)
        return coeff + a + D

    def format_derivatives(self, derivatives: list, variable: str):
        """Given a list of derivatives executes all the derivatives on the variable.

        Parameters
        ----------
        derivatives : list
            list of integers containing the order of the derivative with respect to the
            variable
        variable : str
            variable to be differentiate

        Returns
        -------
        str
            latex code for the derivative.
        """
        var_str = '\\' + variable
        for D, var in zip(derivatives, self.variables):
            for _ in range(D):
                if len(str(variable)) > 1:
                    if '_' in var_str:
                        var_str += "\\" + var
                    else:
                        var_str += '_' + '{' + "\\" + var
                else:
                    if '_' in var_str:
                        var_str += var
                    else:
                        var_str += '_' + '{' + var
        return var_str + '}'

    def format_equation(self, equation: dict) -> str:
        """This code generates the latex string format of the determinant equations

        Parameters
        ----------
        equation : dict
            dictionary containing all the determinant equations

        Returns
        -------
        str
            single string with all the equations typed in latex format
        """
        diff = len(self.constants) - len(self.variables)
        alpha_id = [f'alpha_{idx}' for idx in range(diff)]
        symbolic_constants = self.variables + alpha_id

        latex_equation = ''
        one_term = (len(equation) == 1)
        for term in equation:
            latex_equation += self.dict_to_latex(term, symbolic_constants, one_term) + '+'
        return latex_equation[:-1] + '=0'

    def format_determining_equations(self):
        """Gives the latex code version of remaining determining equation.

        Returns
        -------
        str
            single string with all the equations typed in latex format
        """
        latex_system_of_eqns = ''
        for equation in self.determining_equations.values():
            latex_system_of_eqns += self.format_equation(equation) + "\\" + "\\" + "\n"
        return latex_system_of_eqns

"""Module containing functions for printing in Latex form"""


class Latex():

    def __init__(self, det_eqn, var_dict: dict, var_list: list, constants: list):
        """_summary_

        Parameters
        ----------
        det_eqn : _type_
            _description_
        var_dict : [dict]
            dictionary translating the variable name to latex format
        var_list : [list]
            list with all variables (independent and dependant)list with the constants written in latex format
        constants : [list]
            list containing all constants
        """
        self.det_eqn = det_eqn
        self.var_dict = var_dict
        self.var_list = var_list
        self.constants = constants

    def dict_to_latex(self, term: dict, sym_cte_list: list, one_term: bool):
        """Given a dictionary it returns the symbolic equivalent. It drops all constants if it is  just one term.

        Parameters
        ----------
        term : [dict]
            dictionary containing the information of the term
        sym_cte_list : [list]
            [description]
        one_term : [boolean]
            True if it is the only term in the equation

        Returns
        -------
        [str]
            The latex code to write the term. 
        """
        var = self.var_dict[term['variable']]
        list_devs = term['derivatives']
        cte_power = zip(sym_cte_list, term['constants'])
        a = ''
        if one_term:
            coeff = ''
        else:
            for cte, n in cte_power:
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
        D = self.latex_derivative(list_devs, var)
        latex_term = coeff + a + D
        return latex_term

    def latex_derivative(self, list_devs: list, var: str):
        """Given a list of derivatives executes all the derivatives on the variable.

        Parameters
        ----------
        list_devs : [list]
            list of ints containing the order of the derivative with respect to the variable var_lists.
        var : [str]
            variable to be differentiate

        Returns
        -------
        [str]
            latex code for the derivative.
        """
        D_v = zip(list_devs, self.var_list)
        var_str = '\\' + var
        for D, v in D_v:
            for _ in range(D):
                if len(str(v)) > 1:
                    if '_' in var_str:
                        var_str = var_str + "\\" + v
                    else:
                        var_str = var_str + '_' + '{' + "\\" + v
                else:
                    if '_' in var_str:
                        var_str = var_str + v
                    else:
                        var_str = var_str + '_' + '{' + v
        return var_str + '}'

    def latex_eqn_code(self, eqn: dict) -> str:
        """This code generates the latex string format of the determinant equations

        Parameters
        ----------
        eqn : [dictionary]
            dictionary containing all the determinant equations

        Returns
        -------
        [string]
            single string with all the equations typed in latex format
        """

        sym_cte_list = [var for var in self.var_list]
        for idx in range(len(self.constants) - len(self.var_list)):
            sym_cte_list.append('alpha_' + str(idx))

        A = ''
        for term in eqn:
            one_term = (len(eqn) == 1)
            A += self.dict_to_latex(term, sym_cte_list, one_term) + '+'
        return A[:-1] + '=0'

    def latex_det_eqn(self):
        """Gives the latex code version of remaining determining equation.

        Returns
        -------
        [string]
            single string with all the equations typed in latex format
        """
        latex_code = ''
        for eqn in self.det_eqn.values():
            latex_code += self.latex_eqn_code(eqn) + "\\" + "\\" + "\n"
        return latex_code

"""This file has the structure of the class system. Which has all the information of the physical
system to be analyzed. E.g. the differential equation, rules array, independent and dependent
variables, etc. 
"""

import re
import sympy
from sympy import Derivative as D

from symmetries.utils.combinatorics import list_combinatorics

class System():
    def __init__(self, differential_equation, rules_array:dict,
                 independent_variables:list, dependent_variables:list, constants:list, order:int):
        """Initialization class for system

        Args:
            differential_equation (sympy.add): Differential equation to be analyzed.
            rules_array (dict): Differential equation solved for the higher order derivative.
            independent_variables (list): Independent variables.
            dependent_variables (list): Dependent variables.
            constants (list): Constants of motion of the system.
            order (int): Order of the differential equation.
        """
        self.differential_equation = differential_equation
        self.rules_array = rules_array
        self.independent_variables = independent_variables
        self.dependent_variables = dependent_variables
        self.constants = constants
        self.order = order
        self.n_independent = len(independent_variables)
        self.n_dependent = len(dependent_variables)
        self.infinitesimals: list = []
        self.infinitesimals_dep: list = []
        self.infinitesimals_ind: list = []
        self.dependent_variables_partial_derivatives: list = []
        self.derivatives_subscript_notation: list = []

    def infinitesimals_generator(self):
        """Creates the infinitesimals of the independent and dependents variables. Not the
           derivatives.

        Parameters
        ----------
        independent_variables : list
            list with the independent variables.
        dependent_variables : list
            list with the dependent variables.

        Returns
        -------
        list
            list with the infinitesimals.
        """
        variables = self.independent_variables + self.dependent_variables

        ind_infts = []
        for var in self.independent_variables:
            fun = sympy.Function(f'xi^{var}')
            ind_infts.append(fun(*variables))

        dep_infts = []
        for var in self.dependent_variables:
            fun = sympy.Function(f'eta^{var}'.split('(')[0])
            dep_infts.append(fun(*variables))

        self.infinitesimals = ind_infts+dep_infts
        self.infinitesimals_ind = ind_infts
        self.infinitesimals_dep = dep_infts

    def higher_infinitesimals_generator(self):
        """This functions applies the logic to get the infintesimals of the derivatives.

        Parameters
        ----------
        infinitesimals_of_independent_var : list
            list of the infinitesimals of the independent variables
        infinitesimals_of_dependent_var : list
            list of the infinitesimals of the dependent variables
        order : int
            higher order involved in the system of differential equations
        independent_variables : list
            list with the independent variables
        dependent_variables : list
            list with the dependant variables

        Returns
        -------
        lists
            A list with all possible derivatives of the dependant variables and a list with the 
            respective infinitesimals.
        """
        dep_vars_derivatives = []
        deriv_infints = []
        var_combinatorics = list_combinatorics(self.independent_variables, self.order)
        for deriv_vars_order in var_combinatorics:
            aux_list_deriv = []
            aux_infinitesimals_of_dependent_var = []
            for idx_2, y_aux in enumerate(self.dependent_variables):
                if len(deriv_vars_order) == 1:
                    x_aux =  deriv_vars_order[0]
                    eta_aux = self.infinitesimals_dep[idx_2]
                else:
                    idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                    y_aux = dep_vars_derivatives[idx_1][idx_2]
                    x_aux = deriv_vars_order[-1]
                    eta_aux = deriv_infints[idx_1][idx_2]
                aux_list_deriv.append(D(y_aux, x_aux))
                eta_aux = eta_aux.diff(x_aux)
                for i, ind_i in enumerate(self.independent_variables):
                    eta_aux -= D(y_aux, ind_i)*(
                        self.infinitesimals_ind[i].diff(x_aux))
                aux_infinitesimals_of_dependent_var.append(eta_aux)
            dep_vars_derivatives.append(aux_list_deriv)
            deriv_infints.append(aux_infinitesimals_of_dependent_var)

        infts_dummy = [item for sublist in deriv_infints for item in sublist]
        dep_vars_derivatives = [item for sublist in dep_vars_derivatives for item in sublist]
        self.infinitesimals += infts_dummy
        self.dependent_variables_partial_derivatives = dep_vars_derivatives


    def variable_relabeling(self):
        """Given list of derivatives it changes the partial derivative notation for subscripts
           notation.
        """
        self.derivatives_subscript_notation = []
        for d in self.dependent_variables_partial_derivatives:
            d_str = str(d.args[0]).split('(', maxsplit=1)[0] + '_'
            d_order = list(d.args)
            d_order.pop(0)
            for tup in d_order:
                for _ in range(tup[1]):
                    if '(' in str(tup[0]):
                        v = str(tup[0]).split('(', maxsplit=1)[0]
                    else:
                        v = str(tup[0])
                    d_str = f'{d_str}{v}'
            derivatives_relabel.append(sympy.symbols(d_str))
        self.derivatives_subscript_notation = derivatives_relabel
        self.differential_equation = subs_new_vars(derivatives_relabel,
                                                   self.dependent_variables_partial_derivatives,
                                                   self.differential_equation)


def get_common_factors(XF, dependent_variables:list, independent_variables:list, constants:list):
    """This function creates an empty dictionary where the keys
    are all possible common factors of the determining equations.

    Parameters
    ----------
    XF : symbolic expression
        The Lie operator acting over a differential equation F.
    dependent_variables : list
        list with all dependant variables
    independent_variables : list
        list of all independent variables
    constants : list
        list with all constants

    Returns
    -------
    dict
        empty dictionary where the keys are the possible factorable
        terms for the determining equations.
    """
    S = str(XF.expand())
    S = S.replace(' ', '')
    S = re.sub(r'\*{2}', "&", S)
    S = re.sub(r'Subs\(Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
    S = re.sub(r'Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
    S = re.sub(r'[(eta)(xi)]+\^[\w\(,\)\^]+\**', "", S)
    S = re.sub(r'Derivative', "#", S)
    dep_var_str = [str(ele).replace(' ', '') for ele in dependent_variables]
    for var in dep_var_str:
        var = var.replace("(", "\(").replace(")", "\)")
        S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
    indep_var_str = [str(ele).replace(' ', '') for ele in independent_variables]
    for var in indep_var_str:
        S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
    constants_str = [str(ele).replace(' ', '') for ele in constants]
    for cte in constants_str:
        S = re.sub(f'{cte}\&*\d*', "", S)
    S = re.sub(r'#', "Derivative", S)
    S = re.sub(r'\&', "^", S)
    S = re.sub(r'\*(?=[\+\-])', "", S)
    S = re.sub(r'(?<!\^)\d+\*', '', S)
    S = re.sub(r'(?<=[\+\-])\*+', "", S)
    S = re.sub(r'\*+', "*", S)
    if S[0] == '*':
        S = S[1:]
    keys = re.split('\+|\-', S)
    keys = list(dict.fromkeys(keys))
    if '' in keys:
        keys.remove('')
    return dict.fromkeys(key_ordering(keys))

def key_ordering(keys):
    """Giving a list of strings it organized in a way that the each element is not completely
    inside one of the following terms in the list.

    Parameters
    ----------
    keys : list
        list of no repeated strings

    Returns
    -------
    list
        list with the elements in order
    """
    keys_order = []
    while len(keys_order) < len(keys):
        for key_1 in keys:
            inside = False
            for key_2 in keys:
                if key_2 != key_1 and (key_2 not in keys_order):
                    k1 = key_1.split("*")
                    counter = 0
                    for k in k1:
                        if k in key_2:
                            counter += 1
                    if counter == len(k1):
                        inside = True
                        break
            if not inside and key_1 not in keys_order:
                keys_order.append(key_1)
    return keys_order

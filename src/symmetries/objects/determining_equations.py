"""TODO:Docstring"""
import re
from copy import deepcopy
import sympy
import numpy as np
from symmetries.utils.algebra import key_ordering, str_to_dict, is_zero


class DeterminingEquations():
    """Class"""

    def __init__(self, system, rules_array):
        self.model_info = system
        self.rules_array = rules_array

        self.determining_equations_extended = ""
        self.determining_equations = {}

    def get_group_operator(self):
        """Given a differential equation F, gives the Lie operator acting over F.
        """
        variables = self.model_info.independent_variables \
            + self.model_info.dependent_variables \
            + self.model_info.derivatives_subscript_notation
        l_f = 0
        for var, inft in zip(variables, self.model_info.infinitesimals):
            l_f += inft * \
                sympy.Derivative(self.model_info.differential_equation, var)
        self.determining_equations_extended = sympy.simplify(l_f)

    def variable_relabeling(self):
        """Relabels the determining equation extended version with the derivatives, subscript
           notation

        Parameters
        ----------
        new_labeling : list
            list with the new labeling
        previous_labeling : list
            list with the old symbols
        F : sympy expression
            expression to apply the new labeling

        Returns
        -------
        F: sympy expression
            expression in the new format
        """

        new_labeling = deepcopy(
            self.model_info.dependent_variables_partial_derivatives)
        previous_labeling = deepcopy(
            self.model_info.derivatives_subscript_notation)

        new_labeling.reverse()
        previous_labeling.reverse()

        for new, old in zip(new_labeling, previous_labeling):
            self.determining_equations_extended = \
                self.determining_equations_extended.xreplace({old: new})

    def simplify_rules_array(self):
        """TODO"""
        self.determining_equations_extended = sympy.nsimplify(
            self.determining_equations_extended.subs(self.rules_array))

    def get_common_factors(self):
        """This function adds to the determining equations dictionary the keys of all possible
        common factors of the determining equations.
        """
        S = str(self.determining_equations_extended.expand())
        S = S.replace(' ', '')
        S = re.sub(r'\*{2}', "&", S)
        S = re.sub(r'Subs\(Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
        S = re.sub(r'Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
        S = re.sub(r'[(eta)(xi)]+\^[\w\(,\)\^]+\**', "", S)
        S = re.sub(r'Derivative', "#", S)
        dep_var_str = [str(ele).replace(' ', '')
                       for ele in self.model_info.dependent_variables]
        for var in dep_var_str:
            var = var.replace("(", "\(").replace(")", "\)")
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        indep_var_str = [str(ele).replace(' ', '') for ele in
                         self.model_info.independent_variables]
        for var in indep_var_str:
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        constants_str = [str(ele).replace(' ', '')
                         for ele in self.model_info.constants]
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
        self.determining_equations = dict.fromkeys(key_ordering(keys))

    def construct_determining_equations(self, XF_terms):
        """This function implements the logic to find terms specified
        in dict_det_eqn and group them in that dictionary.

        Parameters
        ----------
        dict_det_eqn : dict
            dictionary where each key has the terms to group by. Each key is a string.
        XF_terms : list
            list of strings where each element is an element of the expanded version an equation.

        Returns
        -------
        dict
            dictionary where in each value contains the factorized terms according to the given
            keys.
        """
        for key in self.determining_equations:
            self.determining_equations[key] = []

        lonely_terms = []
        for element in XF_terms:
            add = False
            for key, value in self.determining_equations.items():
                factors = key.split('*')
                if all(factor in element for factor in factors):
                    for factor in factors:
                        element = element.replace(factor, '')

                    value.append(element)
                    value[-1] = re.sub(r'(?<=[\+\-])\*+', '', value[-1])
                    value[-1] = re.sub(r'(?<=\d|\))\*+', '*', value[-1])
                    value[-1] = re.sub(r'\*\*+', '*', value[-1])

                    if value[-1][-1] == "*":
                        value[-1] = value[-1][:-1]

                    self.determining_equations[key] = value
                    add = True
                    break
            if not add:
                lonely_terms.append(element)
        self.determining_equations['lonely_terms'] = lonely_terms

    def get_determining_equations(self):
        """This function gives the determinant equations in a string format
        """
        S = str(self.determining_equations_extended.expand())
        S = S.replace(' ', '')
        XF_string = re.sub(r'\^', "", S)
        XF_string = re.sub(r'\*{2}', "^", XF_string)
        XF_string = re.sub('\+', ' +', XF_string)
        XF_string = re.sub('\-', ' -', XF_string)
        XF_terms = re.split(' ', XF_string)
        if XF_terms[0] == '':
            XF_terms.pop(0)
        self.construct_determining_equations(XF_terms)

    def encode_determining_equations(self):
        """This function transforms the string version of the determinant equations to a coded
        dictionary format.
        """
        list_var = self.model_info.independent_variables + \
            self.model_info.dependent_variables
        list_all = self.model_info.constants + list_var
        det_eqn = []
        for eqn in self.determining_equations.values():
            aux_list = []
            for str_term in eqn:
                arr_pow = np.zeros(len(list_all))
                arr_deriv = np.zeros(len(list_var))
                term = {"coefficient": 1, "constants": None,
                        "derivatives": None, "variable": None}
                aux_list.append(str_to_dict(sympy.sympify(str_term), term, arr_pow,
                                            arr_deriv, np.array(list_all), np.array(list_var)))
            det_eqn.append(aux_list)
        keys = list(range(len(det_eqn)))
        self.determining_equations = dict(zip(keys, det_eqn))

    def simplify_redundant_equations(self):
        """Given a dict of equations reduces the set by eliminating redundant equations. It is just a
        test to try the logic far from being ready yet
        """
        simplify = True
        exit_param = 0
        zero_terms = {} ## look into this dictionary as to avoid redoing
        det_eqns = self.determining_equations
        det_eqns_aux = deepcopy(self.determining_equations)
        while simplify:
            for idx, eqn in det_eqns.items():
                for _, zero in zero_terms.items():
                    for term in eqn:
                        if is_zero(zero, term) and term in det_eqns_aux[idx]:
                            det_eqns_aux[idx][det_eqns_aux[idx].index(
                                term)] = 0
                            exit_param = 0
                    det_eqns_aux[idx] = list(
                        filter(lambda num: num != 0, det_eqns_aux[idx]))
                    exit_param += 1
                if len(det_eqns_aux[idx]) == 1:
                    zero_terms[idx] = det_eqns_aux[idx][0]
                if exit_param > len(det_eqns):
                    simplify = False
            det_eqns = deepcopy(det_eqns_aux)
        simplify_det_eqns = {k: v for k, v in det_eqns.items() if v}
        for idx in zero_terms:
            simplify_det_eqns[idx] = [zero_terms[idx]]
        return simplify_det_eqns

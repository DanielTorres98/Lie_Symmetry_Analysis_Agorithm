import re
import numpy as np

import sympy
from sympy import Derivative as D
from symmetries.utils.algebra import key_ordering, str_to_dict
from symmetries.utils.symbolic import subs_new_vars

class determining_equations():
    def __init__(self, system):
        self.model_info.system = system
        self.determining_equations_extended = "" 
        self.determining_equations = {}

    def group_operator(self):
        """given a differential equation F, gives the Lie operator acting over F.
        """
        var_inft = zip(self.model_info.independent_variables + self.model_info.dependent_variables
                       + self.model_info.derivatives_subscript_notation,
                       self.model_info.infinitesimals)
        l_f = 0
        for var, inft in var_inft:
            l_f += inft*D(self.model_info.differential_equation, var)
        self.determining_equations_extended = sympy.simplify(l_f)

    def variable_relabeling(self):
        """Relabels the determining equation extended version with the derivatives, subscript
           notation"""
        self.determining_equations_extended = subs_new_vars(
                                        self.model_info.dependent_variables_partial_derivatives,
                                        self.model_info.derivatives_subscript_notation,
                                        self.determining_equations_extended)
    def rules_array_symplification(self):
        self.determining_equations_extended = sympy.nsimplify(
            self.determining_equations_extended.subs(self.model_info.rules_array))

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
        dep_var_str = [str(ele).replace(' ', '') for ele in self.model_info.dependent_variables]
        for var in dep_var_str:
            var = var.replace("(", "\(").replace(")", "\)")
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        indep_var_str = [str(ele).replace(' ', '') for ele in 
                         self.model_info.independent_variables]
        for var in indep_var_str:
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        constants_str = [str(ele).replace(' ', '') for ele in self.model_info.constants]
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
            dictionary where in each value contains the factorized terms according to the given keys.
        """
        for key in self.determining_equations.keys():
            self.determining_equations[key] = []
        lonely_terms = []
        for element in XF_terms:
            add = False
            for key in self.determining_equations:
                factors = key.split('*')
                if len(factors) == 1:
                    if key in element:
                        self.determining_equations[key].append(element.replace(key, ''))
                        self.determining_equations[key][-1] = re.sub(r'(?<=[\+\-])\*+',
                                                    '', self.determining_equations[key][-1])
                        self.determining_equations[key][-1] = re.sub(r'(?<=\d|\))\*+',
                                                    '*', self.determining_equations[key][-1])
                        self.determining_equations[key][-1] = re.sub(r'\*\*+',
                                                    '*', self.determining_equations[key][-1])
                        if self.determining_equations[key][-1][-1] == "*":
                            self.determining_equations[key][-1] = self.determining_equations[
                                                                                        key][-1][:-1]
                        add = True
                        break
                else:
                    inside = True
                    for f in factors:
                        if f not in element:
                            inside = False
                            break
                    if inside:
                        for f in factors:
                            element = element.replace(f, '')
                        self.determining_equations[key].append(element)
                        self.determining_equations[key][-1] = re.sub(r'(?<=[\+\-])\*+',
                                                    '', self.determining_equations[key][-1])
                        self.determining_equations[key][-1] = re.sub(r'(?<=\d)\*+',
                                                    '*', self.determining_equations[key][-1])
                        self.determining_equations[key][-1] = re.sub(r'\*\*+',
                                                    '*', self.determining_equations[key][-1])
                        if self.determining_equations[key][-1][-1] == "*":
                            self.determining_equations[key][-1] = self.determining_equations[
                                                                                        key][-1][:-1]
                        add = True
                        break
            if not add:
                lonely_terms.append(element)
        self.determining_equations['lonely_terms'] = lonely_terms
    
    def get_determining_equations(self):
        """This function gives the determinant equations
        in a string format
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
        self.determining_equations = self.construct_determining_equations(self, XF_terms)

    def encode_determining_equations(self):
        """This function transforms the string version of the determinant equations to a coded 
        dictionary format.
        """
        list_all = (self.model_info.constants + self.model_info.independent_variables + 
                    self.model_info.dependent_variables)
        list_var = self.model_info.independent_variables + self.model_info.dependent_variables
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
        keys = list(np.arange(len(det_eqn)))
        self.determining_equations = dict(zip(keys, det_eqn))


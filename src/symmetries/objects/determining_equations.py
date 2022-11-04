"""TODO:Docstring"""
import re
from copy import deepcopy
import sympy
import numpy as np
from symmetries.utils.algebra import key_ordering, str_to_dict, is_zero
from symmetries.utils.symbolic import get_symbolic_terms, parse_variables
from .system_of_equations import SystemOfEquations
from .system import System


class DeterminingEquations(SystemOfEquations):
    """Class"""

    def __init__(self,
                 rules_array: dict,
                 system: System,
                 independent_variables: list,
                 dependent_variables: list,
                 constants: list,
                 ):

        self.model = system
        self.rules_array = rules_array

        self.determining_equations_extended = ""
        self.determining_equations = {}

        super().__init__(independent_variables, dependent_variables, constants)

        self.general_form = self.obtain_general_form()
        self.parsed_variables = parse_variables(
            independent_variables, dependent_variables)
        self.deleted = {}

    def get_group_operator(self):
        """Given a differential equation F, gives the Lie operator acting over F.
        """
        variables = self.all_variables + self.model.derivatives_subscript_notation
        l_f = 0
        for var, inft in zip(variables, self.model.infinitesimals):
            l_f += inft * \
                sympy.Derivative(self.model.differential_equation, var)
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
            self.model.dependent_variables_partial_derivatives)
        previous_labeling = deepcopy(
            self.model.derivatives_subscript_notation)

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
                       for ele in self.dependent_variables]
        for var in dep_var_str:
            var = var.replace("(", "\(").replace(")", "\)")
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        indep_var_str = [str(ele).replace(' ', '') for ele in
                         self.independent_variables]
        for var in indep_var_str:
            S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
        constants_str = [str(ele).replace(' ', '')
                         for ele in self.constants]
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
        list_var = self.all_variables
        list_all = self.constants + list_var
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

    def obtain_general_form(self):
        """Proposes a general form of solutions for all infinitesimals as a function of
        all the independent and dependent variables accordingly.

        Returns
        -------
        dictionary
            Keys are the infinitesimals eta and eta and values are list of relevant variables.
        """
        general_form = {}
        infinitesimals = self.model.infinitesimals_ind + self.model.infinitesimals_dep

        for inft in infinitesimals:
            l = re.split(r'\W+', str(inft), len(self.all_variables)+1)
            if l[0] not in general_form:
                general_form[l[0]] = {l[1]: deepcopy(self.all_variables)}
            else:
                general_form[l[0]][l[1]] = deepcopy(self.all_variables)

        return general_form

    def simplify_redundant_equations(self):
        """Given a dict of equations reduces the set by eliminating redundant equations. It is just a
        test to try the logic far from being ready yet
        """
        simplify = True
        exit_param = 0
        zero_terms = {}  # look into this dictionary as to avoid redoing
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


    def find_first_derivative_equals_0(self):
        """Finds equations in the determining equations that contain a single term that is a first
        derivative of a given variable and stores it in a dictionary that contains deleted
        dependencies called deleted."""
        for v in self.determining_equations.values():
            if len(v) == 1:
                item = v[0]
                if sum(item['derivatives']) == 1:
                    var = [v for v, order in zip(
                        self.all_variables, item['derivatives']) if order][0]
                    if not ((item['variable'] in self.deleted) and (
                            var in self.deleted[item['variable']])):
                        self.general_form[item['variable'][:-1]
                                          ][item['variable'][-1]].remove(var)
                        if item['variable'] not in self.deleted:
                            self.deleted[item['variable']] = [var]
                            print('deleting', item['variable'], var)
                        else:
                            self.deleted[item['variable']].append(var)
                            print('deleting', item['variable'], var)

    def find_deleted_items_in_equations(self):
        """Searches in the determining equations for equations that contain a term that is an
        already identified first derivative equals 0 or higher order or cross derivative of the
        same."""
        for k, eq in self.determining_equations.items():
            if len(eq) > 1:
                values = deepcopy(eq)
                for item in values:
                    if item['variable'] in self.deleted:
                        variables = [v for v, order in zip(
                            self.all_variables, item['derivatives']) if order]
                        for variable in variables:
                            if (variable in self.deleted[item['variable']]) and (item in eq):
                                print('found', variable, 'in',
                                      item['variable'], 'eq', k)
                                eq.remove(item)

    def delete_second_derivatives(self):
        """Similarly to find_deleted_items_in_equations, the method iterates over the system of
        equations, finds equations containing single higher order derivatives and deletes them
        from the system.
        Note: this is a separate method because it modifies the iterator."""
        det_eqns = deepcopy(self.determining_equations)
        for k, eq in det_eqns.items():
            if len(eq) == 1:
                item = eq[0]  # there is only a single term in the equation
                if item['variable'] in self.deleted:
                    variables = [v for v, order in zip(
                        self.all_variables, item['derivatives']) if order > 1]
                    for variable in variables:
                        if variable in self.deleted[item['variable']]:
                            print('found derivative of', variable,
                                  'in', item['variable'], 'eq', k)
                            del self.determining_equations[k]

    def simplify_iteratively(self):
        """Iterative method, perform the three steps until determining equations does not change.
        """
        while True:
            check_against = deepcopy(self.determining_equations)
            # self.simplify_redundant_equations()
            self.find_first_derivative_equals_0()
            self.find_deleted_items_in_equations()
            self.delete_second_derivatives()

            if check_against == self.determining_equations:
                break

    def print_determining_equations(self):
        """Prints determining equations as symbolic using sympy Matrix."""
        return self._print_symbolic_equations(
            self.determining_equations, 'determining')

    def print_general_form(self):
        """Prints general form as symbolic using sympy Matrix."""
        print('already deleted:', self.deleted)
        print('general form:')
        return self._print_symbolic_equations(
            self.general_form, 'general')

    def _print_symbolic_equations(self, equations, type):
        """Gives the symbolic version of remaining determining equation.

        Parameters
        ----------
        det_eqn : dict
            dictionary will all the determining equations.
        """
        matrix = sympy.Matrix([[]])
        i = 0
        for variable, eqns in equations.items():
            row = None
            if type == 'determining':
                i += 1
                row = sympy.Matrix([[i, sympy.Eq(get_symbolic_terms(
                    eqns, self.parsed_variables, self.constants, self.all_variables), 0)]])
                matrix = matrix.row_insert(i, row)

            elif type == 'general':
                for var, dependencies in eqns.items():
                    row = sympy.Matrix(
                        [f"{variable}({','.join([str(dep) for dep in dependencies])})^{var} "])
                    matrix = matrix.row_insert(i, row)

        return matrix

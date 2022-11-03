
import sympy
import regex as re
from copy import deepcopy
from .system import System
from .determining_equations import DeterminingEquations
from symmetries.utils.symbolic import get_symbolic_terms


class GeneralForm():
    def __init__(self, model: System, system: DeterminingEquations) -> None:
        self.model = model
        self.general_form = self.obtain_general_form()
        self.determining_equations = deepcopy(system.determining_equations)
        self.deleted = {}

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
            l = re.split(r'\W+', str(inft), len(self.model.all_variables)+1)
            if l[0] not in general_form:
                general_form[l[0]] = {l[1]: deepcopy(self.model.all_variables)}
            else:
                general_form[l[0]][l[1]] = deepcopy(self.model.all_variables)

        return general_form

    def find_first_derivative_equals_0(self):
        """Finds equations in the determining equations that contain a single term that is a first
        derivative of a given variable and stores it in a dictionary that contains deleted
        dependencies called deleted."""
        for v in self.determining_equations.values():
            if len(v) == 1:
                l = v[0]
                if sum(l['derivatives']) == 1:
                    var = [self.model.all_variables[n]
                           for n, d in enumerate(l['derivatives']) if d][0]
                    if not (l['variable'] in self.deleted and var in self.deleted[l['variable']]):
                        self.general_form[l['variable'][:-1]
                                          ][l['variable'][-1]].remove(var)
                        if l['variable'] not in self.deleted:
                            self.deleted[l['variable']] = [var]
                            print('deleting', l['variable'], var)
                        else:
                            self.deleted[l['variable']].append(var)
                            print('deleting', l['variable'], var)

    def find_deleted_items_in_equations(self):
        """Searches in the determining equations for equations that contain a term that is an
        already identified first derivative equals 0 or higher order or cross derivative of the
        same."""
        for k, eq in self.determining_equations.items():
            if len(eq) > 1:
                values = deepcopy(eq)
                for item in values:
                    if item['variable'] in self.deleted:
                        for variable, order in zip(self.model.all_variables, item['derivatives']):
                            if (variable in self.deleted[item['variable']]
                            ) and order and (item in eq):
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
                    for variable, order in zip(self.model.all_variables, item['derivatives']):
                        if variable in self.deleted[item['variable']] and order > 1:
                            print('found derivative of', variable,
                                  'in', item['variable'], 'eq', k)
                            del self.determining_equations[k]

    def print_matrix(self):
        print('general form:', self.general_form)  # pretty print this as sympy
        print('already deleted:', self.deleted)
        return self.print_symbolic_equations()

    def print_symbolic_equations(self):
        """Gives the symbolic version of remaining determining equation.

        Parameters
        ----------
        det_eqn : dict
            dictionary will all the determining equations.
        """
        var_dict = self.parse_variables()
        matrix = sympy.Matrix([[]])
        for i, eqn in enumerate(self.determining_equations.values()):
            matrix = matrix.row_insert(i,
                            sympy.Matrix([[i, sympy.Eq(get_symbolic_terms(
                                eqn, var_dict, self.model.constants, self.model.all_variables
                            ), 0)]])
                            )
        return matrix

    def parse_variables(self):
        """Creates dictionary with translation to symbolic terms of variables.

        Returns
        -------
        dictionary
            xit = xi^(t)
        """
        var_dict = {}

        indep_var_str = [str(ele).replace(' ', '')
                         for ele in self.model.independent_variables]
        for dep_var in indep_var_str:
            var = dep_var.split('(')[0]
            if len(var) > 1:
                var_dict[f'xi{var}'] = f'eta^({var})'
            else:
                var_dict[f'xi{var}'] = f'xi^({var})'

        dep_var_str = [str(ele).replace(' ', '')
                       for ele in self.model.dependent_variables]
        for dep_var in dep_var_str:
            var = dep_var.split('(')[0]
            var_dict[f'eta{var}'] = f'eta^({var})'

        return var_dict

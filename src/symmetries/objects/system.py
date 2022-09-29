import sympy
from sympy import Derivative as D

from symmetries.utils.combinatorics import list_combinatorics
from symmetries.utils.symbolic import subs_new_vars

class system():
    def __init__(self, differential_equation, rules_array:dict,
                 independent_variables:list, dependent_variables:list, constants:list, order:int):
        self.differential_equation = differential_equation
        self.rules_array = rules_array
        self.independent_variables = independent_variables
        self.dependent_variables = dependent_variables
        self.constants = constants
        self.order = order
        self.n_independent = len(independent_variables)
        self.n_dependent = len(dependent_variables)
        self.infinitesimals = []
        self.infinitesimals_dep = []
        self.infinitesimals_ind = []
        self.dependent_variables_partial_derivatives = []
        self.derivatives_subscript_notation = []

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
        infts = []
        variables = self.independent_variables + self.dependent_variables
        for var in self.independent_variables:
            fun = sympy.Function(f'xi^{var}')
            infts.append(fun(*variables))
        for var in self.dependent_variables:
            fun = sympy.Function(f'eta^{var}'.split('(')[0])
            infts.append(fun(*variables))
        self.infinitesimals = infts
        self.infinitesimals_of_independent_var = infts[0:self.n_independent]
        self.infinitesimals_of_dependent_var = infts[self.n_independent:(self.n_dependent + 
                                                     self.n_independent)]

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
                    eta_aux = self.infinitesimals_of_dependent_var[idx_2]
                else:
                    idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                    y_aux = dep_vars_derivatives[idx_1][idx_2]
                    x_aux = deriv_vars_order[-1]
                    eta_aux = deriv_infints[idx_1][idx_2]
                aux_list_deriv.append(D(y_aux, x_aux))
                eta_aux = eta_aux.diff(x_aux)
                for i, ind_i in enumerate(self.independent_variables):
                    eta_aux -= D(y_aux, ind_i)*(
                        self.infinitesimals_of_independent_var[i].diff(x_aux))
                aux_infinitesimals_of_dependent_var.append(eta_aux)
            dep_vars_derivatives.append(aux_list_deriv)
            deriv_infints.append(aux_infinitesimals_of_dependent_var)

        infts_dummy = [item for sublist in deriv_infints for item in sublist]
        dep_vars_derivatives = [item for sublist in dep_vars_derivatives for item in sublist]
        self.infinitesimals += infts_dummy
        self.dependent_variables_partial_derivatives = dep_vars_derivatives

    def variable_relabaling(self):
        """Given list of derivatives it changes the partial derivative notation for subscripts
           notation.
        """
        derivatives_relabel = []
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

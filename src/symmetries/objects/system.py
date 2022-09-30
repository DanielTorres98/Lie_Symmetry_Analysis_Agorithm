import sympy
from sympy import Derivative as D

from symmetries.utils.combinatorics import list_combinatorics

class System():
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
            self.derivatives_subscript_notation.append(sympy.symbols(d_str))

    def subs_new_vars(self, new_labeling, previous_labeling):
        """Substitutes notation to another format

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

        new_labeling.reverse()
        previous_labeling.reverse()

        for new, old in zip(new_labeling, previous_labeling):
            self.differential_equation = self.differential_equation.xreplace({old: new})

        new_labeling.reverse()
        previous_labeling.reverse()

    def group_operator(self):
        """given a differential equation F, gives the Lie operator acting over F.

        Parameters
        ----------
        F : sympy expression
            It is the partial differential equation written with sympy symbols
        variables : list
            list containing all variables and the possible derivatives according to the order of the
            differential equation.
        infts : list
            list with the infinitesimals
        """
        all_variables = self.independent_variables \
                        + self.dependent_variables \
                        + self.derivatives_subscript_notation

        l_f = 0
        for var, inft in zip(all_variables, self.infinitesimals):
            l_f += inft*D(self.differential_equation, var)
        self.differential_equation = sympy.simplify(l_f)

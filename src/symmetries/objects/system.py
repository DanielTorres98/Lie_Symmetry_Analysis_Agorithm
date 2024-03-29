"""This file has the structure of the class system. Which has all the information of the physical
system to be analyzed. E.g. the differential equation, rules array, independent and dependent
variables, etc.
"""

from copy import deepcopy
import sympy
from symmetries.utils.combinatorics import list_combinatorics
from .system_of_equations import SystemOfEquations
from typing import Union


class Model(SystemOfEquations):
    """System of equations base class."""

    def __init__(self,
                 differential_equation: Union[sympy.core.add.Add, list],
                 order: Union[int, list],
                 independent_variables: list,
                 dependent_variables: list,
                 constants: list,
                 ):
        """Initialization class for system

        Args:
            differential_equation (sympy.add): Differential equation to be analyzed.
            order (int): Order of the differential equation.
            independent_variables (list): Independent variables.
            dependent_variables (list): Dependent variables.
            constants (list): Constants of motion of the system.
        """
        self.differential_equation = differential_equation
        self.order = max(order) if isinstance(order, list) else order

        self.infinitesimals: list = []
        self.infinitesimals_dep: list = []
        self.infinitesimals_ind: list = []

        self.dependent_variables_partial_derivatives: list = []
        self.derivatives_subscript_notation: list = []

        super().__init__(independent_variables, dependent_variables, constants)

    def infinitesimals_generator(self) -> None:
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
        ind_infts = []
        for var in self.independent_variables:
            fun = sympy.Function(f'xi^{var}')
            ind_infts.append(fun(*self.all_variables))  # pylint: disable=E1102

        dep_infts = []
        for var in self.dependent_variables:
            fun = sympy.Function(f'eta^{var}'.split('(')[0])
            dep_infts.append(fun(*self.all_variables))  # pylint: disable=E1102

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
        var_combinatorics = list_combinatorics(
            self.independent_variables, self.order)
        for deriv_vars_order in var_combinatorics:

            aux_list_deriv = []
            aux_infinitesimals_of_dependent_var = []
            for idx_2, y_aux in enumerate(self.dependent_variables):
                if len(deriv_vars_order) == 1:
                    x_aux = deriv_vars_order[0]
                    # y_aux = y_aux
                    eta_aux = self.infinitesimals_dep[idx_2]

                else:
                    x_aux = deriv_vars_order[-1]
                    idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                    y_aux = dep_vars_derivatives[idx_1][idx_2]
                    eta_aux = deriv_infints[idx_1][idx_2]

                aux_list_deriv.append(sympy.Derivative(y_aux, x_aux))
                eta_aux = eta_aux.diff(x_aux)

                for i, ind_i in enumerate(self.independent_variables):
                    eta_aux -= sympy.Derivative(y_aux, ind_i)*(
                        self.infinitesimals_ind[i].diff(x_aux))

                aux_infinitesimals_of_dependent_var.append(eta_aux)

            dep_vars_derivatives.append(aux_list_deriv)
            deriv_infints.append(aux_infinitesimals_of_dependent_var)

        infts_dummy = [item for sublist in deriv_infints for item in sublist]
        dep_vars_derivatives = [
            item for sublist in dep_vars_derivatives for item in sublist]

        self.infinitesimals += infts_dummy
        self.dependent_variables_partial_derivatives = dep_vars_derivatives

    def variable_relabeling(self):
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
                    v = str(tup[0]).split('(', maxsplit=1)[0]
                    d_str = f'{d_str}{v}'
            derivatives_relabel.append(sympy.symbols(d_str))

        self.derivatives_subscript_notation = derivatives_relabel

        self.differential_equation = self.relabel_differential_eq(
            self.differential_equation)

    def relabel_differential_eq(self, dif_eq: Union[sympy.core.add.Add, list]
                                ) -> Union[sympy.core.add.Add, list]:

        if isinstance(dif_eq, list):
            dif_eq_relabelled = []
            for equation in dif_eq:
                dif_eq_relabelled.append(
                    self.relabel_differential_eq(equation))
            return dif_eq_relabelled

        else:
            assert isinstance(dif_eq, sympy.core.add.Add)
            new_labeling = deepcopy(self.derivatives_subscript_notation)
            previous_labeling = deepcopy(
                self.dependent_variables_partial_derivatives)

            new_labeling.reverse()
            previous_labeling.reverse()

            for new, old in zip(new_labeling, previous_labeling):
                dif_eq = dif_eq.xreplace({old: new})
            for var in self.dependent_variables:
                var_sym = sympy.symbols(str(var).split("(")[0])
                dif_eq = dif_eq.xreplace({var: var_sym})
            return dif_eq

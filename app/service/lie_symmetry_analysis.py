"""This function contains the script to run teh infinitesimals generation, group operator and
simplifying the equations."""
import warnings

from ..src.system import Model
from ..src.determining_equations import DeterminingEquations


def point_symmetries(
    differential_equation,
    order: int,
    f_rules_array: dict,
    independent_variables: list,
    dependent_variables: list,
    constants: list,
):
    # Ignore deprecated warnings. This was implemented to ignore sympy's warning when using a
    # a matrix to show the determining equations.
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    variables = {'independent_variables': independent_variables,
                 'dependent_variables': dependent_variables,
                 'constants': constants}

    model = Model(
        differential_equation=differential_equation,
        order=order,
        **variables
    )
    # This section of the code generates the infinitesimals for all the independent variables and
    # dependent variables.
    #
    model.infinitesimals_generator()
    model.higher_infinitesimals_generator()

    # Relabel derivatives of variables as new symbols to be able to take partial derivatives.
    #
    model.variable_relabeling()

    # Initiating the determining_equations class
    #
    system_of_equations = DeterminingEquations(
        model=model,
        differential_equation=model.differential_equation,
        rules_array=f_rules_array,
        **variables
    )
    system_of_equations.get_common_factors()
    system_of_equations.get_determining_equations()
    system_of_equations.encode_determining_equations()

    system_of_equations.simplify_iteratively()

    return system_of_equations

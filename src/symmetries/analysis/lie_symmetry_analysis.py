"""This function contains the script to run teh infinitesimals generation, group operator and
simplifying the equations."""
import warnings

from symmetries.utils.latex import latex_det_eqn
from symmetries.objects.system import System
from symmetries.objects.determining_equations import DeterminingEquations


def point_symmetries(
    differential_equation,
    order: int,
    f_rules_array: dict,
    independent_variables: list,
    dependent_variables: list,
    constants: list,
    latex: bool = False
):
    # Ignore deprecated warnings. This was implemented to ignore sympy's warning when using a
    # a matrix to show the determining equations.
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    variables = {'independent_variables': independent_variables,
                 'dependent_variables': dependent_variables,
                 'constants': constants}

    model = System(
        differential_equation=differential_equation,
        order=order,
        **variables
    )
    # This section of the code generates the infitesimals for all the independent variables and
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
        system=model,
        rules_array=f_rules_array,
        **variables
    )
    # Applying the group operator over the function F.
    #
    system_of_equations.get_group_operator()
    # Relabel derivatives of variables as new symbols to be able to take partial derivatives.
    
    system_of_equations.variable_relabeling()
    system_of_equations.simplify_rules_array()
    system_of_equations.get_common_factors()
    system_of_equations.get_determining_equations()
    system_of_equations.encode_determining_equations()

    system_of_equations.simplify_iteratively()

    return system_of_equations

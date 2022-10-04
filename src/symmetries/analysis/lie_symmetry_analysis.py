from symmetries.utils.algebra import (get_common_factors, get_det_eqns,
                           simplify_redundant_eqn, str_eqn_to_dict_eqn)
from symmetries.utils.symbolic import (der_relabel, group_operator,
                            higher_infinitesimals_generator,
                            infinitesimals_generator, subs_new_vars,
                            sym_det_eqn)
from symmetries.utils.latex import latex_det_eqn
from symmetries.objects.system import system
from symmetries.objects.determining_equations import determining_equations
import sympy


def point_symmetries(F, order:int, F_rules_array:dict, independent_variables:list,
                     dependent_variables:list, constants:list, latex:bool=False):
    all_variables = independent_variables + dependent_variables
    constants_and_variables = constants + all_variables

    model = system(differential_equation=F, rules_array=F_rules_array,
                   independent_variables=independent_variables,
                   dependent_variables=dependent_variables, constants=constants, order=order)
    # This section of the code generates the infitesimals for all the independent variables and
    # dependent variables.
    #
    model.infinitesimals_generator()
    model.higher_infinitesimals_generator()

    # Relabel derivatives of variables as new symbols to be able to take partial derivatives.
    #
    model.variable_relabaling()

    # Initiating the determining_equations class
    #
    system_of_equations = determining_equations(model)
    # Applying the group operator over the function F.
    #
    system_of_equations.group_operator()
    # Relabel derivatives of variables as new symbols to be able to take partial derivatives.
    #
    system_of_equations.variable_relabeling()
    system_of_equations.rules_array_symplification()
    system_of_equations.get_common_factors()
    system_of_equations.get_determining_equations()
    det_eqn = str_eqn_to_dict_eqn(det_eqn, all_variables, constants_and_variables)
    det_eqn = simplify_redundant_eqn(det_eqn)
    # det_eqns = simplify_redundant_eqn_second_phase(det_eqns)

    # If latex=True prints the latex code for the determining equations.
    #
    # TODO: Fix Latex printing
    #
    if latex:
        backslash_char = "\\"
        latex_dict = {}
        for index, variable in enumerate(independent_variables):
            # If the string length of the variable is bigger than one assumes it is a greek letter.
            #
            if len(str(variable)) > 1:
                latex_dict[f'xse{str(index+1)}'] = f'xi^{"{"}{backslash_char}{variable}'
            else:
                latex_dict[f'xse{str(index+1)}'] = f'xi^{"{"}{variable}{"}"}'
        latex_code = latex_det_eqn(det_eqn, latex_dict, all_variables, constants)
        latex_code = latex_code.replace("+-", "-")
        return print(latex_code)
    return sym_det_eqn(det_eqn, independent_variables, dependent_variables, constants)

"""File with functions regarding symbolic treatment."""
import sympy
from sympy import Derivative as D
from symmetries.utils.combinatorics import list_combinatorics


def dict_to_symb(term, var_dict, var_list, sym_cte_list, one_term):
    """Given a dictionary it returns the symbolic equivalent. It drops all constants if it is
    just one term.

    Parameters
    ----------
    term : dict
        dictionary containing all the information of the term.
    var_dict : dict
        translates the variables from the standard names to users labels.
    var_list : list
        list of strings containing the variables.
    sym_cte_list : list
        list containing the constants in symbolic format.
    one_term : boolean
        True if it is just one term. False otherwise.

    Returns
    -------
    sympy.symbol
        symbolic expression of the derivative.
    """
    var = var_dict[term['variable']]
    list_devs = term['derivatives']
    cte_power = zip(sym_cte_list, term['constants'])
    var_list_str = [str(ele).split('(', maxsplit=1)[0] for ele in var_list]
    a = 1
    if one_term:
        coeff = 1
    else:
        for cte, n in cte_power:
            a *= cte**n
        coeff = term['coefficient']
    derivative = take_derivative(list_devs, var, var_list_str)
    sym_term = coeff*a*derivative
    return sym_term


def take_derivative(dev_list, var_string, var_list):
    """Given a list of derivatives executes all the
       derivatives on the variable.

    Parameters
    ----------
    dev_list : list
        list of ints containing the order of the derivative with respect to the variable var_lists.
    var : symbol
        variable to be differentiate
    var_list : list
        list of independent and dependant
        variables.

    Returns
    -------
    sympy.symbol
        Derivative of the symbolic variable.
    """
    d_var = zip(dev_list, var_list)
    for derivative, var in d_var:
        for _ in range(derivative):
            var_string += '_' + var
    return sympy.symbols(var_string)

def Diff(f, y_var, x_var, order):
    dy = D(y_var, x_var).expand()
    Df = D(f, x_var).expand()
    for _ in range(order):
        Df += dy*D(f, y_var)
        Df = Df.expand()
        y_var = dy
        dy = D(y_var, x_var).expand()
    return Df


def sym_det_eqn(det_eqn, independent_variables, dependent_variables, constants):
    """Gives the symbolic version of remaining
       determining equation.

    Parameters
    ----------
    det_eqn : dict
        dictionary will all the determining equations.
    """
    var_dict = {}
    var_list = independent_variables + dependent_variables

    indep_var_str = [str(ele).replace(' ', '') for ele in independent_variables]
    for dep_var in indep_var_str:
        var = dep_var.split('(')[0]
        if len(var) > 1:
            var_dict[f'xi{var}'] = 'eta^' + \
                '(' + "\\" + var + ')'
        else:
            var_dict[f'xi{var}'] = f'xi^({var})'

    dep_var_str = [str(ele).replace(' ', '') for ele in dependent_variables]
    for dep_var in dep_var_str:
        var = dep_var.split('(')[0]
        if len(var) > 1:
            var_dict[f'eta{var}'] = 'eta^' + \
                '(' + "\\" + var + ')'
        else:
            var_dict[f'eta{var}'] = f'eta^({var})'

    matrix = sympy.Matrix([[]])
    for i, eqn in enumerate(det_eqn.values()):
        matrix = matrix.row_insert(i,
        sympy.Matrix([[i, sympy.Eq(get_symbolic_terms(eqn, var_dict, constants, var_list), 0)]])
        )
    return matrix


def get_symbolic_terms(eqn, var_dict, list_cte, var_list):
    """given a list of dictionaries with the information of each term, returns the symbolic
    equivalent.

    Parameters
    ----------
    eqn : list
        list of dictionaries
    """
    # sym_cte_list = []
    # for idx in range(len(var_list)):
    #     sym_cte_list.append(sympy.symbols(var_list[idx]))
    # for idx in range(len(constants) - len(var_list)):
    #     sym_cte_list.append(sympy.symbols('alpha_' + str(idx)))
    sym_cte_list = list_cte + var_list
    a = 0
    for term in eqn:
        one_term = len(eqn) == 1
        a += dict_to_symb(term, var_dict, var_list,
                             sym_cte_list, one_term)
    return a

def parse_variables(independent_variables, dependent_variables):
    """Creates dictionary with translation to symbolic terms of variables.

    Returns
    -------
    dictionary
        Creates key value pairs passing for example xit to xi^(t).
    """
    var_dict = {}

    indep_var_str = [str(ele).replace(' ', '')
                        for ele in independent_variables]
    for dep_var in indep_var_str:
        var = dep_var.split('(')[0]
        if len(var) > 1:
            var_dict[f'xi{var}'] = f'eta^({var})'
        else:
            var_dict[f'xi{var}'] = f'xi^({var})'

    dep_var_str = [str(ele).replace(' ', '')
                    for ele in dependent_variables]
    for dep_var in dep_var_str:
        var = dep_var.split('(')[0]
        var_dict[f'eta{var}'] = f'eta^({var})'

    return var_dict

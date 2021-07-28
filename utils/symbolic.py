import sympy
from sympy import *

def symbolic_derivative(list_devs, var, var_list):
    """Given a list of derivatives executes all the
       derivatives on the variable.

    Parameters
    ----------
    list_devs : [list]
        list of ints containing the 
        order of the derivative 
        with respect to the variable
        var_lists.
    var : [symbol]
        variable to be differentiate
    var_list : [list]
        list of independant and dependant
        variables.

    Returns
    -------
    [symbol]
        Derivative of the symbolic variable.
    """
    D_v = zip(list_devs, var_list)
    var_str = var
    for D, v in D_v:
        for _ in range(D):
            var_str = var_str + '_' + v
    var = symbols(var_str)
    return var



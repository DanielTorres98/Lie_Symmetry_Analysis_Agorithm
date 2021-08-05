import sympy as sp

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
    var = sp.symbols(var_str)
    return var

def dict_to_symb(term, var_dict, var_list, 
                    sym_cte_list, one_term):
    """Given a dictionary it returns the symbolic
       equivalent. It drops all constants if it is 
       just one term.

    Parameters
    ----------
    term : [dict]
        dictionary containing all thhe information
        of the term.
    var_dict : [dict]
        translates the variables from the standard
        names to users labels.
    var_list : [list]
        list of strings containing the variables.
    sym_cte_list : [list]
        list containing the constants in symbolic
        format.
    one_term : [boolean]
        True if it is just one term. False otherwise.

    Returns
    -------
    [sp.symbol]
        symbolic expression of the derivative.
    """
    var = var_dict[term['variable']]
    list_devs = term['derivatives']
    cte_power = zip(sym_cte_list, term['constants'])
    a = 1
    if one_term:
        coeff = 1
    else:
        for cte, n in cte_power:
            a *= cte**n
        coeff = term['coefficient']
    D = take_derivative(list_devs, var, var_list)
    sym_term = coeff*a*D
    return sym_term

def take_derivative(list_devs, var, var_list):
    """Given a list of derivatives executes all the
       derivatives on the variable.

    Parameters
    ----------
    list_devs : [list]
        list of ints containing the 
        order of the derivative 
        with respect to the variable
        var_lists.
    var : [sp.symbol]
        variable to be differentiate
    var_list : [list]
        list of independant and dependant
        variables.

    Returns
    -------
    [type]
        [description]
    """
    D_v = zip(list_devs, var_list)
    var_str = var
    for D, v in D_v:
        for _ in range(D):
            var_str = var_str + '_' + v
    var = sp.symbols(var_str)
    return var
def dict_to_latex(term, var_dict, var_list, 
                    sym_cte_list, one_term):
    """Given a dictionary it returns the symbolic
       equivalent. It drops all constants if it is 
       just one term.

    Parameters
    ----------
    term : [dict]
        dictionary containing the information of the term
    var_dict : [dict]
        dictionary translating the variable name to 
        latex format
    var_list : list with all variables (independant and 
        dependant)
        list with the constants writen in latex format
    sym_cte_list : [list]
        [description]
    one_term : [boolean]
        True if it is the only term in the equation

    Returns
    -------
    [str]
        The latex code to write the term. 
    """
    var = var_dict[term['variable']]
    list_devs = term['derivatives']
    cte_power = zip(sym_cte_list, term['constants'])
    a = ''
    if one_term:
        coeff = ''
    else:
        for cte, n in cte_power:
            if len(cte)>1:
                if n == 1:
                    a = a + "\\" + str(cte)
                if n > 1:
                    a = a + "\\" + str(cte) + '^' + str(n) 
            else:
                if n == 1:
                    a = a + str(cte)
                if n > 1:
                    a = a + str(cte) + '^' + str(n)           
        coeff = str(term['coefficient'])
    D = latex_derivative(list_devs, var, var_list)
    latex_term = coeff + a + D
    return latex_term

def latex_derivative(list_devs, var, var_list):
    """Given a list of derivatives executes all the
       derivatives on the variable.

    Parameters
    ----------
    list_devs : [list]
        list of ints containing the 
        order of the derivative 
        with respect to the variable
        var_lists.
    var : [str]
        variable to be differentiate
    var_list : [list]
        list of independant and dependant
                          variables.

    Returns
    -------
    [str]
        latex code for the derivative.
    """
    D_v = zip(list_devs, var_list)
    var_str = '\\' + var
    for D, v in D_v:
        for _ in range(D):
            if len(v)>1:
                if '_' in var_str:
                    var_str = var_str + "\\" +  v
                else:
                    var_str = var_str + '_' + '{' + "\\" + v
            else:
                if '_' in var_str:
                    var_str = var_str + v
                else:
                    var_str = var_str + '_' + '{' + v
    var = var_str + '}'
    return var
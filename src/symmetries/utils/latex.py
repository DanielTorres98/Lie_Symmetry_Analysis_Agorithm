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

def latex_eqn_code(eqn, var_dict, var_list, constants):
    """This code generates the latex string format of 
    the determinant equations

    Parameters
    ----------
    eqn : [dictionary]
        dictionary containing all the determinant eqautions
    var_dict : [dictionary]
        dictionary mapping the name of the variable to it's 
        latex code 
    var_list : [list]
        list with all independant and dependant variables in
        string format
    constants : [list]
        list containing all constants

    Returns
    -------
    [string]
        returns a single string with  all the equations 
        typed in latex format
    """
    sym_cte_list = []
    for idx in range(len(var_list)):
        sym_cte_list.append(var_list[idx])
    for idx in range(len(constants) - len(var_list)):
        sym_cte_list.append('alpha_' + str(idx))
    A = ''
    for term in eqn:
        one_term = False
        if len(eqn) == 1:
            one_term = True
        A = A + dict_to_latex(term, var_dict, var_list,
                                 sym_cte_list, one_term) + '+'
    A = A[:-1]
    A = A + '=0'
    return A


def latex_det_eqn(det_eqn, var_dict, var_list, constants):
    """Gives the latex code version of remaining
       determining equation.

       Args:
       det_eqn (dict): dictionary will all the
                       determining equations. 
    """
    latex_code = ''
    for eqn in det_eqn.values():
        latex_code = latex_code + latex_eqn_code(eqn,
                                                 var_dict, var_list, constants) + "\\" + "\\" + "\n"
    return latex_code
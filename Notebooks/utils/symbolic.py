import sympy as sp
from sympy import Derivative as D
from utils.combinatorics import list_combinatorics


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
    [sp.symbol]
        symbolic representation of the derivative
    """
    D_v = zip(list_devs, var_list)
    var_str = var
    for D, v in D_v:
        for _ in range(D):
            var_str = var_str + '_' + v
    var = sp.symbols(var_str)
    return var


def infinitesimals_generator(list_indep, list_dep):
    """Creates the infinitesimals of the independent and
    depentants variables. Not the derivatives. 

    Parameters
    ----------
    list_indep : [list]
        list with the independent variables.
    list_dep : [list]
        list with the dependent variables.

    Returns
    -------
    [list]
        list with the infinitesimals. 
    """
    infts = []
    variables = list_indep + list_dep
    for var in list_indep:
        f = f'xi^{var}'
        infts.append(sp.Function(f)(*variables))
    for var in list_dep:
        f = f'eta^{var}'.split('(')[0]
        infts.append(sp.Function(f)(*variables))
    return infts


def higher_infinitesimals_generator(list_inft_indep, list_inft,
                                    order, list_indep, list_dep):
    """This funtions applies the logic to get the infintesimals
       of the derivatives.

    Parameters
    ----------
    list_inft_indep : [list]
        list of the infinitesimals of the independent variables
    list_inft : [list]
        list of the infinitesimals of the dependent variables
    order : [int]
        higher order involved in the system of differential equations
    list_indep : [list]
        list with the independant variables
    list_dep : [type]
        list with the dependant variables

    Returns
    -------
    [lists]
        A list with all possible derivatives of the dependant
        variables and a list with the respective infinitesimals
    """
    dep_vars_derivatives = []
    inft_derivatives = []
    var_combinatorics = list_combinatorics(list_indep, order)
    for deriv_vars_order in var_combinatorics:
        aux_list_deriv = []
        aux_list_inft = []
        for idx_2 in range(len(list_dep)):
            if len(deriv_vars_order) == 1:
                y_aux = list_dep[idx_2]
                x_aux =  deriv_vars_order[0]
                dummy_eta = list_inft[idx_2]
            else:
                idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                y_aux = dep_vars_derivatives[idx_1][idx_2]
                x_aux = deriv_vars_order[-1]
                dummy_eta = inft_derivatives[idx_1][idx_2]
            aux_list_deriv.append(D(y_aux, x_aux))
            dummy_eta = dummy_eta.diff(x_aux)           
            for i in range(len(list_indep)):
                dummy_eta -= D(y_aux, list_indep[i])*(
                    list_inft_indep[i].diff(x_aux))
            aux_list_inft.append(dummy_eta)
        dep_vars_derivatives.append(aux_list_deriv)
        inft_derivatives.append(aux_list_inft)
    return inft_derivatives, dep_vars_derivatives

def group_operator(F, variables, infts):
    """given a differential equation F, gives the
    Lie operator acting over F. 

    Parameters
    ----------
    F : [sympy expression]
        It is the partial differential equation written with
        sympy symbols
    variables : [list]
        list containing all variables and the possible derivatives
        according to the order of the differential equation.
    infts : [list]
        list with the infinitesimals
    """
    var_inft = zip(variables, infts)
    LF = 0
    for var, inft in var_inft:
        LF += sp.simplify(inft*D(F, var))
    return LF

def der_relabel(dep_vars_derivatives, F):
    derivatives_relabel = []
    for d in dep_vars_derivatives:
        d_str = str(d.args[0]).split('(')[0] + '_' 
        d_order = list(d.args)
        d_order.pop(0)
        for tup in d_order:
            for t in range(tup[1]):
                d_str = f'{d_str}{tup[0]}'
        derivatives_relabel.append(sp.symbols(d_str))

    derivatives_relabel.reverse()
    dep_vars_derivatives.reverse()
    names_derivatives = zip(derivatives_relabel, 
                            dep_vars_derivatives)

    for name, deriv in names_derivatives:
        F = F.subs({deriv:name})
    derivatives_relabel.reverse()
    dep_vars_derivatives.reverse()
    return F, derivatives_relabel
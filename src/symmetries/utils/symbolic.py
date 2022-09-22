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
    var_list_str = [str(ele).split('(')[0] for ele in var_list]
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


def take_derivative(dev_list, var, var_list):
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
            var_string = var_string + '_' + var
    return sympy.symbols(var_string)


def infinitesimals_generator(list_indep, list_dep):
    """Creates the infinitesimals of the independent and dependents variables. Not the derivatives.

    Parameters
    ----------
    list_indep : list
        list with the independent variables.
    list_dep : list
        list with the dependent variables.

    Returns
    -------
    list
        list with the infinitesimals.
    """
    infts = []
    variables = list_indep + list_dep
    for var in list_indep:
        fun = sympy.Function(f'xi^{var}')
        infts.append(fun(*variables))
    for var in list_dep:
        fun = sympy.Function(f'eta^{var}'.split('(')[0])
        infts.append(fun(*variables))
    return infts


def higher_infinitesimals_generator(list_inft_indep, list_inft,
                                    order, list_indep, list_dep):
    """This functions applies the logic to get the infintesimals of the derivatives.

    Parameters
    ----------
    list_inft_indep : list
        list of the infinitesimals of the independent variables
    list_inft : list
        list of the infinitesimals of the dependent variables
    order : int
        higher order involved in the system of differential equations
    list_indep : list
        list with the independent variables
    list_dep : list
        list with the dependant variables

    Returns
    -------
    lists
        A list with all possible derivatives of the dependant variables and a list with the 
        respective infinitesimals.
    """
    dep_vars_derivatives = []
    deriv_infints = []
    var_combinatorics = list_combinatorics(list_indep, order)
    for deriv_vars_order in var_combinatorics:
        aux_list_deriv = []
        aux_list_inft = []
        for idx_2, y_aux in enumerate(list_dep):
            if len(deriv_vars_order) == 1:
                x_aux =  deriv_vars_order[0]
                eta_aux = list_inft[idx_2]
            else:
                idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                y_aux = dep_vars_derivatives[idx_1][idx_2]
                x_aux = deriv_vars_order[-1]
                eta_aux = deriv_infints[idx_1][idx_2]
            aux_list_deriv.append(D(y_aux, x_aux))
            eta_aux = eta_aux.diff(x_aux)
            for i, ind_i in enumerate(list_indep):
                eta_aux -= D(y_aux, ind_i)*(
                    list_inft_indep[i].diff(x_aux))
            aux_list_inft.append(eta_aux)
        dep_vars_derivatives.append(aux_list_deriv)
        deriv_infints.append(aux_list_inft)
    return deriv_infints, dep_vars_derivatives


def group_operator(F, variables, infts):
    """given a differential equation F, gives the Lie operator acting over F.

    Parameters
    ----------
    F : sympy expression
        It is the partial differential equation written with sympy symbols
    variables : list
        list containing all variables and the possible derivatives according to the order of the
        differential equation.
    infts : list
        list with the infinitesimals
    """
    var_inft = zip(variables, infts)
    l_f = 0
    for var, inft in var_inft:
        l_f += inft*D(F, var)
    return sympy.simplify(l_f)

def der_relabel(dep_vars_derivatives, F):
    """Given list of derivatives it changes the partial derivative notation for subscripts.

    Parameters
    ----------
    dep_vars_derivatives : list
        list of all possible derivatives
    F : sympy expression
        Functional to substitute in.

    Returns
    -------
    F : sympy expression
        The functional with the new notation of the derivatives.

    derivatives_relabel : list
        A list with the new symbols for the visualization.
    """
    derivatives_relabel = []
    for d in dep_vars_derivatives:
        d_str = str(d.args[0]).split('(')[0] + '_'
        d_order = list(d.args)
        d_order.pop(0)
        for tup in d_order:
            for _ in range(tup[1]):
                if '(' in str(tup[0]):
                    v = str(tup[0]).split('(')[0]
                else:
                    v = str(tup[0])
                d_str = f'{d_str}{v}'
        derivatives_relabel.append(sympy.symbols(d_str))

    F = subs_new_vars(derivatives_relabel, dep_vars_derivatives, F)
    return F, derivatives_relabel

def subs_new_vars(new_labeling, previous_labeling, F):
    """Substitutes notation to another format

    Parameters
    ----------
    new_labeling : list
        list with the new labeling
    previous_labeling : list
        list with the old symbols
    F : sympy expression
        expression to apply the new labeling

    Returns
    -------
    F: sympy expression
        expression in the new format
    """
    new_labeling.reverse()
    previous_labeling.reverse()

    for new, old in zip(new_labeling, previous_labeling):
        F = F.xreplace({old: new})
    new_labeling.reverse()
    previous_labeling.reverse()
    return F

def deriv_infts(infts, variables, order):
    """_summary_

    Parameters
    ----------
    infts : _type_
        _description_
    variables : _type_
        _description_
    order : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    deriv_infitnits = []
    for deriv_vars_order in list_combinatorics(variables, order):
        for inft in infts:
            for var in deriv_vars_order:
                inft_aux = D(inft, var)
            deriv_infitnits.append(inft_aux)
    return deriv_infitnits

def Diff(f, y_var, x_var, order):
    dy = D(y_var, x_var).expand()
    Df = D(f, x_var).expand()
    for _ in range(order):
        Df += dy*D(f, y_var)
        Df = Df.expand()
        y_var = dy
        dy = D(y_var, x_var).expand()
    return Df

def higher_infinitesimals_generator_2(list_inft_indep, list_inft,
                                    order, list_indep, list_dep):
    """This functions applies the logic to get the infintesimals
       of the derivatives.

    Parameters
    ----------
    list_inft_indep : list
        list of the infinitesimals of the independent variables
    list_inft : list
        list of the infinitesimals of the dependent variables
    order : int
        higher order involved in the system of differential equations
    list_indep : list
        list with the independent variables
    list_dep : type
        list with the dependant variables

    Returns
    -------
    lists
        A list with all possible derivatives of the dependant
        variables and a list with the respective infinitesimals.
    """
    dep_vars_derivatives = []
    deriv_infints = []
    var_combinatorics = list_combinatorics(list_indep, order)
    for deriv_vars_order in var_combinatorics:
        aux_list_deriv = []
        aux_list_inft = []
        for idx_2, dep in enumerate(list_dep):
            if len(deriv_vars_order) == 1:
                y_aux = dep
                x_aux =  deriv_vars_order[0]
                eta_aux = list_inft[idx_2]
            else:
                idx_1 = var_combinatorics.index(deriv_vars_order[:-1])
                y_aux = dep_vars_derivatives[idx_1][idx_2]
                x_aux = deriv_vars_order[-1]
                eta_aux = deriv_infints[idx_1][idx_2]
            aux_list_deriv.append(D(y_aux, x_aux))
            eta_aux = eta_aux.diff(x_aux)

            for i, ind in enumerate(list_indep):
                eta_aux -= D(y_aux, ind)*(list_inft_indep[i].diff(x_aux))
            aux_list_inft.append(eta_aux)
        dep_vars_derivatives.append(aux_list_deriv)
        deriv_infints.append(aux_list_inft)

    return deriv_infints, dep_vars_derivatives

def sym_det_eqn(det_eqn, list_indep, list_dep, constants):
    """Gives the symbolic version of remaining
       determining equation.

    Parameters
    ----------
    det_eqn : dict
        dictionary will all the determining equations.
    """
    var_dict = {}
    var_list = list_indep + list_dep

    indep_var_str = [str(ele).replace(' ', '') for ele in list_indep]
    for dep_var in indep_var_str:
        var = dep_var.split('(')[0]
        if len(var) > 1:
            var_dict[f'xi{var}'] = 'eta^' + \
                '(' + "\\" + var + ')'
        else:
            var_dict[f'xi{var}'] = f'xi^({var})'

    dep_var_str = [str(ele).replace(' ', '') for ele in list_dep]
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
        one_term = False
        if len(eqn) == 1:
            one_term = True
        a += dict_to_symb(term, var_dict, var_list,
                             sym_cte_list, one_term)
    return a

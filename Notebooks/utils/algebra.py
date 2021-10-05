import copy
import re

import numpy as np
import sympy as sp


def is_zero(zero_term, term):
    """Given a term that is zero, returns true if
       term is zero as well.

    Parameters
    ----------
    zero_term : [dict]
        a dictionary containing the
        information of the zero term
    term : [dict]
        a dictionary containing the
        information of the term

    Returns
    -------
    [boolean]
        Returns True if the term is zero as well, False
        otherwise.
    """
    return zero_term['variable'] == term['variable'] and\
        compare_derivatives(zero_term['derivatives'], term['derivatives'])


def compare_derivatives(D1, D2):
    """Given to lists with the information of the derivatives
       tells if the second term contains a derivative equal 
       or higher order for all possible derivatives.

    Parameters
    ----------
    D1 : [list]
        list of the order of derivatives
        for each variable for term 1.
    D2 : [list]
        list of the order of derivatives
        for each variable for term 1.

    Returns
    -------
    [Boolean]
        Returns False if at least one derivative in D2 is
        of a lower order than in D1. Returns True otherwise.
    """
    D1D2 = zip(D1, D2)
    for d1, d2 in D1D2:
        if d2 < d1:
            return False
    return True


def drop_constants(eqn):
    """If a term has the same constants it drops them.

    Parameters
    ----------
    eqn : [list]
        list containing all terms

    Returns
    -------
    [list]
        A list containing all terms but without the constants
        multiplying the whole equation. 
    """
    equal_terms = True
    equal_constants = True
    terms_info = eqn[0]["constants"]
    coeff_info = [abs(eqn[0]["coefficient"]), eqn[0]["constants"]]
    for term in eqn:
        if coeff_info != \
                [abs(term["coefficient"]), term["constants"]]:
            equal_terms = False
        if terms_info != term["constants"]:
            equal_constants = False
    if equal_constants:
        N = len(eqn[0]["constants"])
        for i in range(len(eqn)):
            if equal_terms:
                eqn[i]["coefficient"] = int(eqn[i]["coefficient"]
                                            / abs(eqn[i]["coefficient"]))
            eqn[i]["constants"] = [0 for i in range(N)]
    return eqn


def get_common_factors(XF, list_dep, list_indep, constants):
    """This funtion creates an empty dictionary where the keys
    are all possible commmon factors of the determining equations.

    Parameters
    ----------
    XF : [symbolic expresion]
        The Lie operator acting over a differential equation F.
    list_dep : [list]
        list with all dependant variables
    list_indep : [list]
        list of all independant varables
    constants : [list]
        list with all constants

    Returns
    -------
    [dict]
        empty dictionary where the keys are the possible factorizable
        terms for the determining equations.
    """
    S = str(XF.expand())
    S = S.replace(' ', '')
    S = re.sub(r'\*{2}', "&", S)
    S = re.sub(r'Subs\(Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
    S = re.sub(r'Derivative\([(eta)(xi)\^][\w\(,\)\^]*\&*\d*', "", S)
    S = re.sub(r'[(eta)(xi)]+\^[\w\(,\)\^]+\**', "", S)
    S = re.sub(r'Derivative', "#", S)
    dep_var_str = [str(ele).replace(' ', '') for ele in list_dep]
    for var in dep_var_str:
        var = var.replace("(", "\(").replace(")", "\)")
        S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
    indep_var_str = [str(ele).replace(' ', '') for ele in list_indep]
    for var in indep_var_str:
        S = re.sub(f'(?<!\(|,){var}\&*\d*', "", S)
    constants_str = [str(ele).replace(' ', '') for ele in constants]
    for cte in constants_str:
        S = re.sub(f'{cte}\&*\d*', "", S)
    S = re.sub(r'#', "Derivative", S)
    S = re.sub(r'\&', "^", S)
    S = re.sub(r'\*(?=[\+\-])', "", S)
    S = re.sub(r'(?<!\^)\d+\*', '', S)
    S = re.sub(r'(?<=[\+\-])\*+', "", S)
    S = re.sub(r'\*+', "*", S)
    keys = re.split('\+|\-', S)
    keys = list(dict.fromkeys(keys))
    if '' in keys:
        keys.remove('')
    return dict.fromkeys(key_ordering(keys))


def key_ordering(keys):
    """Giving a list of strings it organized in a way
    that the each element is not completely inside 
    one of the following terms in the list.

    Parameters
    ----------
    keys : [list]
        list of no repeated strings

    Returns
    -------
    [list]
        list with the elements in order
    """
    keys_order = []
    while len(keys_order) < len(keys):
        for key_1 in keys:
            inside = False
            for key_2 in keys:
                if key_2 != key_1 and (key_2 not in keys_order):
                    k1 = key_1.split("*")
                    counter = 0
                    for k in k1:
                        if k in key_2:
                            counter += 1
                    if counter == len(k1):
                        inside = True
                        break
            if not inside and key_1 not in keys_order:
                keys_order.append(key_1)
    return keys_order


def get_det_eqns(XF, dict_det_eqn):
    """This function gives the determinant equations
    in a string format

    Parameters
    ----------
    XF : [sp.add]
        symbolic representation of the Lie operator acting 
        over the differential equation. At this step XF already
        used the rules array. 
    dict_det_eqn : [dict]
        dictionary with all the derivatives terms in XF

    Returns
    -------
    [dict]
        dictionary with the string version of the determinant 
        equations. Each key has the information of the 
        term that was grouped by and equated to zero. 
    """
    S = str(XF.expand())
    S = S.replace(' ', '')
    XF_string = re.sub(r'\^', "", S)
    XF_string = re.sub(r'\*{2}', "^", XF_string)
    XF_string = re.sub('\+', ' +', XF_string)
    XF_string = re.sub('\-', ' -', XF_string)
    XF_terms = re.split(' ', XF_string)
    if XF_terms[0] == '':
        XF_terms.pop(0)
    return group_terms(dict_det_eqn, XF_terms)


def group_terms(dict_det_eqn, XF_terms):
    """This function implements the logic to find terms specified
    in dict_det_eqn and group them in that dictionary.

    Parameters
    ----------
    dict_det_eqn : [dict]
        dictionary where each key has the terms to group by. Each 
        key is a string.
    XF_terms : [list]
        list of strings where each element is an element of the 
        expanded version an equation.

    Returns
    -------
    [dict]
        dictionary where in each value contains the factorized terms
        according to the given keys. 
    """
    for key in dict_det_eqn.keys():
        dict_det_eqn[key] = []
    lonely_terms = []
    for element in XF_terms:
        add = False
        for key in dict_det_eqn:
            factors = key.split('*')
            if len(factors) == 1:
                if key in element:
                    dict_det_eqn[key].append(element.replace(key, ''))
                    dict_det_eqn[key][-1] = re.sub(r'(?<=[\+\-])\*+',
                                                   '', dict_det_eqn[key][-1])
                    dict_det_eqn[key][-1] = re.sub(r'(?<=\d|\))\*+',
                                                   '*', dict_det_eqn[key][-1])
                    dict_det_eqn[key][-1] = re.sub(r'\*\*+',
                                                   '*', dict_det_eqn[key][-1])
                    if dict_det_eqn[key][-1][-1] == "*":
                        dict_det_eqn[key][-1] = dict_det_eqn[key][-1][:-1]
                    add = True
                    break
            else:
                inside = True
                for f in factors:
                    if f not in element:
                        inside = False
                        break
                if inside:
                    for f in factors:
                        element = element.replace(f, '')
                    dict_det_eqn[key].append(element)
                    dict_det_eqn[key][-1] = re.sub(r'(?<=[\+\-])\*+',
                                                   '', dict_det_eqn[key][-1])
                    dict_det_eqn[key][-1] = re.sub(r'(?<=\d)\*+',
                                                   '*', dict_det_eqn[key][-1])
                    dict_det_eqn[key][-1] = re.sub(r'\*\*+',
                                                   '*', dict_det_eqn[key][-1])
                    if dict_det_eqn[key][-1][-1] == "*":
                        dict_det_eqn[key][-1] = dict_det_eqn[key][-1][:-1]
                    add = True
                    break
        if not add:
            lonely_terms.append(element)
    dict_det_eqn['lonely_terms'] = lonely_terms
    return dict_det_eqn


def str_eqn_to_dict_eqn(dict_det_eqn, list_var, list_all):
    """This function transforms the string version of the 
    determinant equations to dictionary format.

    Parameters
    ----------
    dict_det_eqn : [dict]
        dictionary containing all the determinant equations
    list_var : [list]
        list with all the variables
    list_all : [list]
        list with all the variables and constants (constants
        go first).

    Returns
    -------
    [dict]
        a dictionary of lists where each element of the list is 
        another dictionary with the information of each term 
        of the determining equations of the system. 
    """
    det_eqn = []
    for eqn in dict_det_eqn.values():
        aux_list = []
        print(eqn)
        for str_term in eqn:
            print(str_term)
            arr_pow = np.zeros(len(list_all))
            arr_deriv = np.zeros(len(list_var))
            term = {"coefficient": 1, "constants": None,
                    "derivatives": None, "variable": None}
            aux_list.append(str_to_dict(sp.sympify(str_term), term, arr_pow,
                  arr_deriv, np.array(list_all), np.array(list_var)))
        det_eqn.append(aux_list)
    keys = list(np.arange(len(det_eqn)))
    return dict(zip(keys, det_eqn))


def str_to_dict(f, term, arr_pow, arr_deriv, list_all, list_var):
    """Returns a dictionary that saves all the information of a 
    symbolic expresion in a determining equation.

    Parameters
    ----------
    f : [simpy expression]
        A simbolic expresion to analyze
    term : [dict]
        dictionary containing all the information of the term
    arr_pow : [numpy array]
        array with the power of each constant or variable
        multiplying the term.
    arr_deriv : [numpy array]
        array with the order of the derivative with respect to
         the variables in list_var
    list_all : [list]
        list with all constants and variables. Constants go first
    list_var : [list]
        list with all dependant and independant variables. 

    Returns
    -------
    [dict]
        dictionary with the information of the symbolic term in 
        a determining equation. 
    """
    if f.is_Mul:
        for i in f.args:
            str_to_dict(i, term, arr_pow, arr_deriv, list_all, list_var)
    else:
        if f.args == ():
            if f.is_Integer:
                term['coefficient'] = f
            else:
                idx = np.where(list_all == f)
                p = 1
                arr_pow[idx] = p

        if f.is_Function:
            idx = np.where(list_all == f)
            p = 1
            arr_pow[idx] = p
            if idx[0].size == 0:
                term['variable'] = str(f).split('(')[0]

        if f.is_Pow:
            var = f.args[0]
            if var.is_Derivative:
                term['variable'] = var.args[0]
                for v in var.args[1:]:
                    idx = np.where(list_var == v[0])
                    arr_deriv[idx] = v[1]
            else:
                idx = np.where(list_all == f.args[0])
                arr_pow[idx] = f.args[1]

        if type(f) == type(sp.Subs(list_var[0],
                                   list_var[0], list_var[0])):
            s = str(f)
            if s.endswith('))'):
                if '(' in re.split(',', s)[-1]:
                    subs = re.split(',',s)[-2].strip(' ')
                    subs_t = re.split(',',s)[-1].strip(')').strip(' ')    # Changed
                    subs_t = subs_t + ')'
                else:
                    subs = re.split(',', s)[-3].strip(' ')
                    subs_t = re.split(',', s)[-2] + ',' + re.split(',', s)[-1]
                    subs_t = subs_t[:-1].strip(' ')
            else:
                subs = re.split(',',s)[-2].strip(' ')
                subs_t = re.split(',',s)[-1].strip(')').strip(' ')
            if ')' in subs:
                subs = subs.strip(')')
            s = s.replace(subs, subs_t)
            f = sp.sympify(parens(s).strip('('))
            if isinstance(f, tuple):
                f = f[0]

        if f.is_Derivative:
            term['variable'] = str(f.args[0]).split("(")[0]
            for v in f.args[1:]:
                idx = np.where(list_var == v[0])
                arr_deriv[idx] = v[1]

    term['constants'] = list(arr_pow.astype(int))
    term['derivatives'] = list(arr_deriv.astype(int))
    return term


def parens(s):
    i = s.count(')') - 1
    groups = s[s.find('('):].split(')')
    return (')'.join(groups[:i]) + ')')

def simplify_redundant_eqn(det_eqns):
    """given a dict of equations reduces the set 
       by eliminating redundant eqns. It is just a test
       to try the logic far from being ready yet

       Args:
       det_eqns (dict): 
    """
    simplify = True
    N = len(det_eqns)
    exit = 0
    zero_terms = {}
    det_eqns_aux = copy.deepcopy(det_eqns)
    while simplify:
        for idx, eqn in det_eqns.items():
            for i, zero in zero_terms.items():
                for term in eqn:
                    if is_zero(zero, term) and term in det_eqns_aux[idx]:
                        det_eqns_aux[idx][det_eqns_aux[idx].index(term)] = 0
                        exit = 0
                det_eqns_aux[idx] = list(
                    filter(lambda num: num != 0, det_eqns_aux[idx]))
                exit += 1
            if len(det_eqns_aux[idx]) == 1:
                zero_terms[idx] = det_eqns_aux[idx][0]
            if exit > N:
                simplify = False
        det_eqns = copy.deepcopy(det_eqns_aux)
    simplify_det_eqns = {k: v for k, v in det_eqns.items() if v}
    for idx in zero_terms:
        simplify_det_eqns[idx] = [zero_terms[idx]]
    return simplify_det_eqns

def simplify_redundant_eqn_second_phase(det_eqn):
    for idx, eqn in det_eqn.items():
        det_eqn[idx] = de.drop_constants(eqn)
    return det_eqn

def drop_constants(eqn):
    """If a term has the same constants it drops them.

    Parameters
    ----------
    eqn : [list]
        list containing all terms

    Returns
    -------
    [list]
        A list containing all terms but without the constants
        multiplying the whole equation. 
    """
    equal_terms = True
    equal_constants = True
    terms_info = eqn[0]["constants"]
    coeff_info = [abs(eqn[0]["coefficient"]), eqn[0]["constants"]]
    for term in eqn:
        if coeff_info != \
                [abs(term["coefficient"]), term["constants"]]:
            equal_terms = False
        if terms_info != term["constants"]:
            equal_constants = False
            return eqn
    if equal_constants:
        N = len(eqn[0]["constants"])
        for i in range(len(eqn)):
            if equal_terms:
                eqn[i]["coefficient"] = int(eqn[i]["coefficient"]
                                            / abs(eqn[i]["coefficient"]))
            eqn[i]["constants"] = [0 for i in range(N)]
    return eqn

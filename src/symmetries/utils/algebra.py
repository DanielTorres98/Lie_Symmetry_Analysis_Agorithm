"""This module contains functions regarding the algebra."""
import re

import numpy as np
import sympy

from symmetries.utils.DEint import compare_derivatives, drop_constants

def is_zero(zero_term:dict, term:dict):
    """Given a term that is zero, returns true if
       term is zero as well.

    Parameters
    ----------
    zero_term : dict
        a dictionary containing the
        information of the zero term
    term : dict
        a dictionary containing the
        information of the term

    Returns
    -------
    boolean
        Returns True if the term is zero as well, False
        otherwise.
    """
    return zero_term['variable'] == term['variable'] and\
        compare_derivatives(zero_term['derivatives'], term['derivatives'])



def key_ordering(keys):
    """Giving a list of strings it organized in a way that the each element is not completely
    inside one of the following terms in the list.

    Parameters
    ----------
    keys : list
        list of no repeated strings

    Returns
    -------
    list
        list with the elements in order
    """
    keys_order = []
    while len(keys_order) < len(keys):
        filtered_keys = [key for key in keys if key not in keys_order]
        for key_1 in filtered_keys:
            inside = False
            for key_2 in filtered_keys:
                if key_2 != key_1 and all(k in key_2 for k in key_1.split("*")):
                    inside = True
                    break
            if not inside:
                keys_order.append(key_1)
    return keys_order



def str_eqn_to_dict_eqn(dict_det_eqn, list_var, list_all):
    """This function transforms the string version of the determinant equations to dictionary
    format.

    Parameters
    ----------
    dict_det_eqn : dict
        dictionary containing all the determinant equations
    list_var : list
        list with all the variables
    list_all : list
        list with all the variables and constants (constants
        go first).

    Returns
    -------
    dict
        a dictionary of lists where each element of the list is another dictionary with the
        information of each term of the determining equations of the system.
    """
    det_eqn = []
    for eqn in dict_det_eqn.values():
        aux_list = []
        for str_term in eqn:
            arr_pow = np.zeros(len(list_all))
            arr_deriv = np.zeros(len(list_var))
            term = {"coefficient": 1, "constants": None,
                    "derivatives": None, "variable": None}
            aux_list.append(str_to_dict(sympy.sympify(str_term), term, arr_pow,
                  arr_deriv, np.array(list_all), np.array(list_var)))
        det_eqn.append(aux_list)
    keys = list(np.arange(len(det_eqn)))
    return dict(zip(keys, det_eqn))


def str_to_dict(f, term, arr_pow, arr_deriv, list_all, list_var):
    """Returns a dictionary that saves all the information of a symbolic expression in a
    determining equation.

    Parameters
    ----------
    f : sympy expression
        A symbolic expression to analyze
    term : dict
        dictionary containing all the information of the term
    arr_pow : numpy array
        array with the power of each constant or variable
        multiplying the term.
    arr_deriv : numpy array
        array with the order of the derivative with respect to
         the variables in list_var
    list_all : list
        list with all constants and variables. Constants go first
    list_var : list
        list with all dependant and independent variables.

    Returns
    -------
    dict
        dictionary with the information of the symbolic term in a determining equation.
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
                for var in var.args[1:]:
                    idx = np.where(list_var == var[0])
                    arr_deriv[idx] = var[1]
            else:
                idx = np.where(list_all == f.args[0])
                arr_pow[idx] = f.args[1]

        if type(f) == type(sympy.Subs(list_var[0], list_var[0], list_var[0])):
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
            f = sympy.sympify(parens(s).strip('('))
            if isinstance(f, tuple):
                f = f[0]

        if f.is_Derivative:
            term['variable'] = str(f.args[0]).split("(")[0]
            for var in f.args[1:]:
                idx = np.where(list_var == var[0])
                arr_deriv[idx] = var[1]

    term['constants'] = list(arr_pow.astype(int))
    term['derivatives'] = list(arr_deriv.astype(int))
    return term

def parens(s):
    i = s.count(')') - 1
    groups = s[s.find('('):].split(')')
    return ')'.join(groups[:i]) + ')'

def simplify_redundant_eqn_second_phase(det_eqn):
    for idx, eqn in det_eqn.items():
        det_eqn[idx] = drop_constants(eqn)
    return det_eqn
    

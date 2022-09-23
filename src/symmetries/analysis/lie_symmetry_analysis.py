from symmetries.utils.algebra import (get_common_factors, get_det_eqns,
                           simplify_redundant_eqn, str_eqn_to_dict_eqn)
from symmetries.utils.symbolic import (der_relabel, group_operator,
                            higher_infinitesimals_generator,
                            infinitesimals_generator, subs_new_vars,
                            sym_det_eqn)
from symmetries.utils.latex import latex_det_eqn
import sympy as sp


def point_symmetries(F, order, F_rules_array, list_indep, list_dep, list_cte, latex=False):
    list_var = list_indep + list_dep
    list_all = list_cte + list_indep + list_dep
    latex_dict = {}
    # for i in range(0, N_indep):
    #     if len(var_list[i]) > 1:
    #         latex_dict['xse' + str(i+1)] = 'xi^' + \
    #             '{' + "\\" + var_list[i] + '}'
    #     else:
    #         latex_dict['xse' + str(i+1)] = 'xi^' + "{" + var_list[i] + '}'
    infts = infinitesimals_generator(list_indep, list_dep)
    n_indep = len(list_indep)
    n_dep = len(list_dep)
    infts_ind = infts[0:n_indep]
    infts_dep = infts[n_indep:(n_dep + n_indep)]
    inft_derivatives, dep_vars_derivatives = higher_infinitesimals_generator(infts_ind, infts_dep,
                                                                             order, list_indep, list_dep)
    infts_dummy = [item for sublist in inft_derivatives for item in sublist]
    dep_vars_derivatives = [item for sublist in dep_vars_derivatives for item in sublist]
    infts = infts + infts_dummy
    F, deriv_names = der_relabel(dep_vars_derivatives, F)
    vars_and_derivatives = list_indep + list_dep + deriv_names
    XF = group_operator(F, vars_and_derivatives, infts)
    XF = subs_new_vars(dep_vars_derivatives, deriv_names, XF)
    XF = sp.simplify(XF.subs(F_rules_array))
    empty_det_eqn = get_common_factors(XF, list_dep, list_indep, list_cte)
    det_eqn = get_det_eqns(XF, empty_det_eqn)
    det_eqn = str_eqn_to_dict_eqn(det_eqn, list_var, list_all)
    det_eqn = simplify_redundant_eqn(det_eqn)
    # det_eqns = simplify_redundant_eqn_second_phase(det_eqns)
    # if latex:
    #     latex_code = latex_det_eqn(det_eqns,
    #                     latex_dict, list_var, list_cte)
    #     latex_code = latex_code.replace("+-", "-")
    #     return print(latex_code)
    return sym_det_eqn(det_eqn, list_indep, list_dep, list_cte)

# from utils.algebra import simplify_redundant_eqn
# from utils.symbolic import sym_det_eqn

# def point_symmetries_det_eqn(cvs_path, N_indep, N_dep,
#                              var_list, cte_list, simplify=True, latex=False):
#     """Gives a symbolic representation of all equations
#        from a cvs file.

#        Args:
#        cvs_path (str): path where the cvs is.
#        N_indep (int): number of independent variables.
#        N_dep (int): number of dependant variables.
#        var_list (list): list of all variables labels
#                         in string format.
#        simplify (boolean): If True, simplifies the set
#                            of determining equations.
#     """
#     var_dict = {}
#     for i in range(0, N_indep):
#         var_dict['xse' + str(i+1)] = 'xi^' + '(' + var_list[i] + ')'
#     for i in range(0, N_dep):
#         k = i + N_indep
#         if len(var_list[k]) > 1:
#             var_dict['eta' + str(i+1)] = 'eta^' + \
#                 '(' + "\\" + var_list[k] + ')'
#         else:
#             var_dict['eta' + str(i+1)] = 'eta^' + "(" + var_list[k] + ')'
#     latex_dict = {}
#     for i in range(0, N_indep):
#         if len(var_list[i]) > 1:
#             latex_dict['xse' + str(i+1)] = 'xi^' + \
#                 '{' + "\\" + var_list[i] + '}'
#         else:
#             latex_dict['xse' + str(i+1)] = 'xi^' + "{" + var_list[i] + '}'
#     for i in range(0, N_dep):
#         k = i + N_indep
#         if len(var_list[k]) > 1:
#             latex_dict['eta' + str(i+1)] = 'eta^' + \
#                 '{' + "\\" + var_list[k] + '}'
#         else:
#             latex_dict['eta' + str(i+1)] = 'eta^' + "{" + var_list[k] + '}'
#     det_eqns, cte_list = cvs_to_list(cvs_path, var_list)
#     if simplify:
#         det_eqns = simplify_redundant_eqn(det_eqns)
#     #   det_eqns = simplify_redundant_eqn_second_phase(det_eqns)
#     # if latex:
#     #     latex_code = latex_det_eqn(det_eqns,
#     #                                latex_dict, var_list, cte_list)
#     #     latex_code = latex_code.replace("+-", "-")
#     #     return print(latex_code)
#     return sym_det_eqn(det_eqns, var_dict, var_list, cte_list), det_eqns

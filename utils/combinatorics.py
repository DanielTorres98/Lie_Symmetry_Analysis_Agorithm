import itertools

def list_combinatorics(variables):
    """Given a list of variables returns
       all possible combinations.

    Parameters
    ----------
    variables : [list]
        A list with the symbolic variables.

    Returns
    -------
    [list]
        A list of lists with all possible 
        combinations
    """
    var_combinatorics = []
    for L in range(0, len(variables)+1):
        for subset in itertools.combinations(variables, L):
            var_combinatorics.append(list(subset))
    return var_combinatorics

def higher_infinitesimals_generator(list_inft_indep ,list_inft, Order, list_indep, list_dep):
    dep_vars_derivatives = []
    inft_derivatives = []
    var_combinatorics = cb.list_combinatorics(list_indep)
    depvars_etas = zip(list_dep, list_inft)
    indvar_xis =  zip(list_indep, list_inft_indep)
    for ind_vars in var_combinatorics:
        for dep_var, eta in depvars_etas: 
            if len(ind_vars) == 1:
                y = dep_var
                x = ind_vars[0]
                e = eta
            else:
                idx = var_combinatorics.index(ind_vars[:-1])  
                y = dep_vars_derivatives[idx] 
                x = ind_vars[-1]
                e = inft_derivatives[idx]
            dep_vars_derivatives.append(D(y, x))
            dummy = e.diff(x)
            for ind_var, xi in indvar_xis:
                dummy += D(y, ind_var)*(xi.diff(x))
            inft_derivatives.append(dummy)
    return dep_vars_derivatives, inft_derivatives
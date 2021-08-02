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
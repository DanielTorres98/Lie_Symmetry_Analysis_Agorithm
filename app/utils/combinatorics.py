"""list_combinatorics function"""

from itertools import combinations_with_replacement as cwr

def list_combinatorics(variables, order):
    """Given a list of variables returns all possible combinations.

    Parameters
    ----------
    variables : list
        A list with the symbolic variables.

    order : int
        order of the differential equation

    Returns
    -------
    list
        A list of lists with all possible combinations
    """
    var_combinatorics = []
    for order_i in range(1, order+1):
        var_combinatorics = var_combinatorics + list(cwr(variables, order_i))
    var_combinatorics = [list(x) for x in var_combinatorics]
    return var_combinatorics

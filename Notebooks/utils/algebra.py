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
    return zero_term['variable']==term['variable'] and\
          compare_derivatives(zero_term['derivatives'], term['derivatives'])
    
def compare_derivatives(D1,D2):
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
    D1D2 =  zip(D1, D2)
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
        if  coeff_info != \
            [abs(term["coefficient"]), term["constants"]]:
            equal_terms = False
        if terms_info != term["constants"]:
            equal_constants = False
    if equal_constants:
        N = len(eqn[0]["constants"])
        for i in range(len(eqn)):
            if equal_terms:
                eqn[i]["coefficient"] = int(eqn[i]["coefficient"]
                /abs(eqn[i]["coefficient"]))
            eqn[i]["constants"] = [0 for i in range(N)]
    return eqn

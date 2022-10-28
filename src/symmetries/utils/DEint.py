"""I want to put the interpretation code I write in a separate file,
 so I can call them when I need to in a main analysis file"""
import pandas as pd

def int_list(string_list:str):
    """ Receives a list written in a string format and transforms it into a list object.

    Parameters
    ----------
    s : str
        string version of a list.
    """
    string_list = string_list.strip('[')
    string_list = string_list.strip(']')
    interpreted = [int(n) for n in string_list.split(', ')]
    return interpreted

def get_terms(equation: str):
    """A function that will split the equations into terms. The terms are always separated
    by a + or -.

    Parameters
    ----------
    equation : str
        equation written in a string form.
    """
    equation = equation.replace('+','-')
    
    # Break the equation into it's terms
    #
    term_list = equation.split('-')
    # Now strip the whitespace
    #
    terms = [term.strip(' ') for term in term_list]
    if terms[0] == '':
        terms.pop(0)
    # This gets rid of the empty string that results from the equation starting with
    #  a negative sign.
    return terms


def get_signs(equation: str):
    """For a given equation in a string format will give a list with the sign of each term.

    Parameters
    ----------
    equation : str
        equation written in a string form.
    """
    sign_array = []
    if equation[0] != "-":
        sign_array.append(1)
    for char in equation:
        if char == "+":
            sign_array.append(1)
        elif char == '-':
            sign_array.append(-1)
    return sign_array


def find_constants(data: pd.DataFrame, var_list: list):
    """Find the constants in the equations in the data set.

    Parameters
    ----------
    data : pd.DataFrame
        pd.DataFrame containing all equations.
    """
    var_csv_label_list = []
    for k in range(1, len(var_list)+1):
        var_csv_label_list.append('z' + str(k))

    for equation in data:
        for term in get_terms(equation):
            parts = term.split('*')
            if len(parts) > 1:
                try:
                    parts[0] = int(parts[0])
                except:
                # Here we catch any greek letters, and then put an integer of 1 at the beginning
                # of the list
                    parts[0] = parts[0].strip('(')
                    parts = [1] + parts

            if len(parts) > 2:
                # At this point we definitely have a greek letter: print(parts[1:-1])
                for letter in parts[1:-1]:
                    if '^' in letter:
                        # If the constant is raised to a power
                        nc = letter.split('^')[0]
                    else:
                        nc = letter
                    if nc not in var_csv_label_list:
                        var_csv_label_list.append(nc)
    return var_csv_label_list

def process_term(term, constant_list, var_list):
    """A function that splits a single term up into individual parts

    Parameters
    ----------
    term : str
        Individual term of the equation in string format
    constant_list : list
        list with the constants of the set of original differential equations
    """
    # First need to separate coefficients
    #
    parts = term.split('*')

    # Make the numerical coefficient into an integer
    #
    if len(parts) > 1:
        try:
            parts[0] = int(parts[0])
        except:
            # Here we catch any greek letters, and then put an integer of 1 at the beginning
            # of the list
            #
            parts[0] = parts[0].strip('(')
            parts = [1] + parts

    # Need to strip in 2 stages as to not delete the function information
    #
    var_str = ''
    for k in range(1, len(var_list)):
        var_str = var_str + 'z' + str(k) + ', '
    var_str = var_str + 'z' + str(k + 1) + '])'
    parts[-1] = parts[-1].strip(var_str).strip('[')
    greek_list = []
    if len(parts) > 2:
        # At this point we definitely have a greek letter: print(parts[1:-1])
        #
        for constant in constant_list:
            i = 0
            # In general we have something like 'η^2'.
            # We need to add the power of each greek letter to the list.
            for letter in parts[1:-1]:
                if constant in letter:
                    if '^' in letter:
                        i += int(letter.split('^')[1])
                    else:
                        i += 1
                else:
                    continue
            greek_list.append(i)
    else:
        greek_list = [0]*len(constant_list)

    fnc, derivatives = process_derivative(parts[-1], var_list)
    return [parts[0], greek_list, derivatives, fnc]

# Function to process the derivatives.

def process_derivative(derivative, var_list):
    """ Takes the derivative and variable information from a string

    Parameters
    ----------
    derivative : str
        A containing the information of which variable is being derived and the order
        of the derivatives
    """
    if 'Derivative' in derivative:
        derivative = derivative.strip('Derivative')
        ints = derivative[:derivative.find(']')+1]
        fnc = derivative[derivative.find(']')+1:].strip('[]')
    else:
        ints = [0 for i in var_list]
        return derivative, ints
    return fnc, ints(ints)

def term_to_dict(term):
    """ Takes the list with the information of the each term and put it in a dictionary

    Parameters
    ----------
    term : list
        List containing the information of an individual term.
    """
    if not isinstance(term[0], int):
        term[0] = 1
    term_dict = {"coefficient": term[0],
                 "constants": term[1],
                 "derivatives": term[2],
                 "variable": term[3]}
    return term_dict

def eqn_process(equation, constant_list, var_list):
    """ Takes an equation and split it in a list of dictionaries
        saving all the relevant information for each term

    Parameters
    ----------
    term : list
        List containing the information of an individual term.
    """
    signs = get_signs(equation)
    terms = get_terms(equation)
    signs_terms =  zip(signs, terms)
    list_terms = []
    for sign, term in signs_terms:
        processed_term = process_term(term, constant_list, var_list)
        term_dict = term_to_dict(processed_term)
        term_dict["coefficient"] = term_dict["coefficient"]*sign
        list_terms.append(term_dict)
    return list_terms

def is_zero(zero_term, term):
    """Given a term that is zero, returns true if term is zero as well.

    Parameters
    ----------
    zero_term: dict
        a dictionary containing the information of the zero term
    term : dict
        a dictionary containing the information of the term
    """
    return zero_term['variable']==term['variable'] and\
          compare_derivatives(zero_term['derivatives'], term['derivatives'])


def compare_derivatives(der_1:list, der_2:list)-> bool:
    """Given to lists with the information of the derivatives tells if the second term contains a
     derivative equal or higher order for all possible derivatives.

    Parameters
    ----------
    der_1 : list
        list of the order of derivatives
        for each variable for term 1.
    der_2 : list
        list of the order of derivatives
        for each variable for term 1.

    Returns
    -------
    Boolean
        Returns False if at least one derivative in D2 is
        of a lower order than in D1. Returns True otherwise.
    """
    for d_1, d_2 in zip(der_1, der_2):
        if d_2 < d_1:
            return False
    return True

def drop_constants(eqn: list[dict]):
    """If a term has the same constants it drops them.

    Parameters
    ----------
    eqn : list
        dict containing all terms

    Returns
    -------
    list
        A list containing all terms but without the constants multiplying the whole equation.
    """
    equal_terms = True
    equal_constants = True
    terms_info = eqn[0]["constants"]
    coefficient_info = [abs(eqn[0]["coefficient"]), eqn[0]["constants"]]

    for term in eqn:
        if coefficient_info != [abs(term["coefficient"]), term["constants"]]:
            equal_terms = False
        if terms_info != term["constants"]:
            equal_constants = False

    if equal_constants:
        for i, _ in enumerate(eqn):
            if equal_terms:
                eqn[i]["coefficient"] = int(eqn[i]["coefficient"]/ abs(eqn[i]["coefficient"]))
            eqn[i]["constants"] = [0 ]*len(eqn[0]["constants"])
    return eqn

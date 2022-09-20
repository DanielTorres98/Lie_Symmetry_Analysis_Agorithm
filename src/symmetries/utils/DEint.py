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
    term_list = equation.split('-')
    # Now strip the whitespace
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
        else:
            continue
    return sign_array


def find_constants(data: pd.DataFrame, var_list: list):
    """Find the constants in the equations in the data set.

    Parameters
    ----------
    data : pd.DataFrame
        dataframe containing all equations.
    """
    var_csv_label_list = []
    for k in range(1, len(var_list)):
        var_csv_label_list.append('z' + str(k))
    var_csv_label_list.append('z' + str(k+1))
    cnst_list = var_csv_label_list
    for equation in data:
        for term in get_terms(equation):
            parts = term.split('*')
            if len(parts) > 1:
                try:
                    parts[0] = int(parts[0])
                except:
                # Here we catch any greek letters, and then put an integer of 1 at the beginning of the list
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
                    if nc in cnst_list:
                        continue
                    else:
                        cnst_list.append(nc)
    return cnst_list
                                                

def process_term(term, cnst_list, var_list):
    """ # A function that splits a single term up into individual parts

        Parameters
        term (str): Individual term of the equation in
                    string format
        cnst_list (list): list with the constants of the
                          set of original differential equations
    """
    # First need to separate coefficients
    parts = term.split('*')
    # Make the numerical coefficient into an integer
    if len(parts) > 1:
        try:
            parts[0] = int(parts[0])
        except:
        # Here we catch any greek letters, and then put an integer of 1 at the beginning of the list
            parts[0] = parts[0].strip('(')
            parts = [1] + parts

    # Need to strip in 2 stages as to not delete the function information
    var_str = ''
    for k in range(1, len(var_list)):
        var_str = var_str + 'z' + str(k) + ', '
    var_str = var_str + 'z' + str(k + 1) + '])'
    parts[-1] = parts[-1].strip(var_str).strip('[')
    greek_list = []
    if len(parts) > 2:
        # At this point we definitely have a greek letter: print(parts[1:-1])
        for cnst in cnst_list:
            i = 0
            # In general we have something like 'η^2'. 
            # We need to add the power of each greek letter to the list.
            for letter in parts[1:-1]:
                if cnst in letter:
                    if '^' in letter:
                        i += int(letter.split('^')[1])
                    else:
                        i += 1
                else:
                    continue
            greek_list.append(i)          
    else:
        greek_list = [0 for cnst in cnst_list]
        pass
    fnc, derivatives = process_derivative(parts[-1], var_list)
    return [parts[0], greek_list, derivatives, fnc]

# Function to process the derivatives.

def process_derivative(derivative, var_list):
    """ Takes the derivative and variable information from a string

        Parameters
        derivative (str): A containing the information of which
                          variable is being derivend and the order
                          of the derivatives
    """
    if 'Derivative' in derivative:
        derivative = derivative.strip('Derivative')
        int_list = derivative[:derivative.find(']')+1]
        fnc = derivative[derivative.find(']')+1:].strip('[]')
    else:
        int_list = [0 for i in var_list]
        return derivative, int_list
    return fnc, int_list(int_list)

def term_to_dict(term):
    """ Takes the list with the information of the each term and put it in a dictionary 

        Parameters
        term (list): List containing the information of
                     an individual term. 
    """
    if type(term[0]) != int:
        term[0] = 1
    term_dict = {"coefficient": term[0],
                 "constants": term[1],
                 "derivatives": term[2],
                 "variable": term[3]}
    return term_dict

def eqn_process(equation, cnst_list, var_list):
    """ Takes an equation and split it in a list of dictionaries
        saving all the relevant information for each term

        Parameters
        term (list): List containing the information of
                     an individual term. 
    """
    signs = get_signs(equation)
    terms = get_terms(equation)
    signs_terms =  zip(signs, terms)
    list_terms = []
    for sign, term in signs_terms:
        processed_term = process_term(term, cnst_list, var_list)
        term_dict = term_to_dict(processed_term)
        term_dict["coefficient"] = term_dict["coefficient"]*sign
        list_terms.append(term_dict)
    return list_terms

def is_zero(zero_term, term):
    """Given a term that is zero, returns true if term is zero as well.

        Parameters
        zero_term (dict): a dictionary containing the
                          information of the zero term
        term (dict):      a dictionary containing the
                          information of the term
    """
    return zero_term['variable']==term['variable'] and\
          compare_derivatives(zero_term['derivatives'], term['derivatives'])


def compare_derivatives(D1,D2):
    """Given to lists with the information of the derivatives
       tells if the second term contains a derivative equal 
       or higher order for all possible derivatives.

        Parameters
        D1 (list): list of derivatives of term 1
        D2 (list): list of derivatives of term 1
    """
    D1D2 =  zip(D1, D2)
    for d1, d2 in D1D2:
        if d2 < d1:
            return False
    return True

def drop_constants(eqn):
    equal_terms = True
    coeff_info = [abs(eqn[0]["coefficient"]), eqn[0]["constants"]]
    for term in eqn:
        if  coeff_info != \
            [abs(term["coefficient"]), term["constants"]]:
            equal_terms = False
    if equal_terms:
        N = len(eqn[0]["constants"])
        for i in range(len(eqn)):
            eqn[i]["coefficient"] = int(eqn[i]["coefficient"]
            /abs(eqn[i]["coefficient"]))
            eqn[i]["constants"] = [0 for i in range(N)]
    return eqn
def dict_to_symb(term, var_dict, var_list, 
                    sym_cte_list, one_term):
    """Given a dictionary it returns the symbolic equivalent. It drops all constants if it is
       just one term.

        Parameters
        eqn (dict): dictionary with the information of
                    the term.
    """
    var = var_dict[term['variable']]
    list_devs = term['derivatives']
    cte_power = zip(sym_cte_list, term['constants'])
    a = 1
    if one_term:
        coeff = 1
    else:
        for cte, n in cte_power:
            a *= cte**n
        coeff = term['coefficient']
    D = take_derivative(list_devs, var, var_list)
    sym_term = coeff*a*D
    return sym_term

def dict_to_latex(term, var_dict, var_list, 
                    sym_cte_list, one_term):
    """Given a dictionary it returns the symbolic equivalent. It drops all constants if it is
       just one term.

        Parameters
        eqn (dict): dictionary with the information of
                    the term.
    """
    var = var_dict[term['variable']]
    list_devs = term['derivatives']
    cte_power = zip(sym_cte_list, term['constants'])
    a = ''
    if one_term:
        coeff = ''
    else:
        for cte, n in cte_power:
            if len(cte)>1:
                if n == 1:
                    a = a + "\\" + str(cte)
                if n > 1:
                    a = a + "\\" + str(cte) + '^' + str(n) 
            else:
                if n == 1:
                    a = a + str(cte)
                if n > 1:
                    a = a + str(cte) + '^' + str(n)           
        coeff = str(term['coefficient'])
    D = latex_derivative(list_devs, var, var_list)
    latex_term = coeff + a + D
    return latex_term

def take_derivative(list_devs, var, var_list):
    """Given a list of derivatives executes all the derivatives on the variable.

        Parameters
        list_devs (list): list of ints containing the 
                          order of the derivative 
                          with respect to the variable
                          var_lists.
        var (symbol):     variable to be differentiate
        var_list (list):  list of independant and dependant
                          variables.
    """
    D_v = zip(list_devs, var_list)
    var_str = var
    for D, v in D_v:
        for _ in range(D):
            var_str = var_str + '_' + v
    var = symbols(var_str)
    return var

def latex_derivative(list_devs, var, var_list):
    """Given a list of derivatives executes all the derivatives on the variable.

        Parameters
        list_devs (list): list of ints containing the 
                          order of the derivative 
                          with respect to the variable
                          var_lists.
        var (symbol):     variable to be differentiate
        var_list (list):  list of independant and dependant
                          variables.
    """
    D_v = zip(list_devs, var_list)
    var_str = '\\' + var
    for D, v in D_v:
        for _ in range(D):
            if len(v)>1:
                if '_' in var_str:
                    var_str = var_str + "\\" +  v
                else:
                    var_str = var_str + '_' + '{' + "\\" + v
            else:
                if '_' in var_str:
                    var_str = var_str + v
                else:
                    var_str = var_str + '_' + '{' + v
    var = var_str + '}'
    return var
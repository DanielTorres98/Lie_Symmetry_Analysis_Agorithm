from symmetries.analysis.lie_symmetry_analysis import point_symmetries
import sympy as sp
from sympy import Derivative as D


def test_heat_equation():
    x, t, k = sp.symbols('x t k')
    u = sp.Function('u')(x, t)
    independent_variables = [x, t]
    dependent_variables = [u]
    constants = [k]
    F = D(u, t) - k*D(u, x, x)
    F_rules_array = {D(u, x, x): 1/k*D(u, t)}
    order = 2

    system_of_equations = point_symmetries(
        F, order, F_rules_array, independent_variables, dependent_variables, constants)

    assert str(system_of_equations.general_form
        ) == "{'xi': {'x': [x, t], 't': [t]}, 'eta': {'u': [x, t, u(x, t)]}}"


def test_schrödinger_equation():
    x, t, k, i = sp.symbols('x t k i')
    u = sp.Function('u')(x, t)
    independent_variables = [x, t]
    dependent_variables = [u]
    constants = [k, i]
    F = -D(u, x, x) + x**2*u-i*k*D(u, t)
    F_rules_array = {D(u, x, x): x**2-i*k*D(u, t)}
    order = 2

    system_of_equations = point_symmetries(
        F, order, F_rules_array, independent_variables, dependent_variables, constants)

    assert str(system_of_equations.general_form
        ) == "{'xi': {'x': [x, t], 't': [t]}, 'eta': {'u': [x, t, u(x, t)]}}"

def test_blasius_equation():
    x = sp.symbols('x')
    y = sp.Function('y')(x)
    independent_variables = [x]
    dependent_variables = [y]
    constants = []
    F =  D(y,x,x,x) + y*D(y,x,x)
    F_rules_array = {D(y,x,x,x):-y*D(y,x,x)}
    order = 3

    system_of_equations = point_symmetries(
        F, order, F_rules_array, independent_variables, dependent_variables, constants)

    assert str(system_of_equations.general_form
        ) == "{'xi': {'x': [x]}, 'eta': {'y': [x, y(x)]}}"

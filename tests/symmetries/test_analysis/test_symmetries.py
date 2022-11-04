from symmetries.analysis.lie_symmetry_analysis import point_symmetries
import sympy as sp
from sympy import Derivative as D
from symmetries.objects.general_form import GeneralForm


def test_heat_equation():
    x, t, k = sp.symbols('x t k')
    u = sp.Function('u')(x, t)
    independent_variables = [x, t]
    dependent_variables = [u]
    constants = [k]
    F = D(u, t) - k*D(u, x, x)
    F_rules_array = {D(u, x, x): 1/k*D(u, t)}
    order = 2
    F.expand()

    model, system_of_equations = point_symmetries(
        F, order, F_rules_array, independent_variables, dependent_variables, constants)
    general_form = GeneralForm(model, system_of_equations)
    general_form.simplify_iteratively()
    assert str(general_form.general_form) == "{'xi': {'x': [x, t], 't': [t]}, 'eta': {'u': [x, t, u(x, t)]}}"

def test_schrödinger_equation():
    x, t, k, i = sp.symbols('x t k i')
    u = sp.Function('u')(x, t)
    list_indep = [x, t]
    list_dep = [u]
    list_cte = [k, i]
    F =  -D(u,x,x) +x**2*u-i*k*D(u,t)
    F_rules_array = {D(u,x,x):x**2-i*k*D(u,t)}
    order = 2

    model, system_of_equations = point_symmetries(F, order, F_rules_array, list_indep, list_dep, list_cte)

    general_form = GeneralForm(model, system_of_equations)
    general_form.simplify_iteratively()

    assert str(general_form.general_form) == "{'xi': {'x': [x, t], 't': [t]}, 'eta': {'u': [x, t, u(x, t)]}}"
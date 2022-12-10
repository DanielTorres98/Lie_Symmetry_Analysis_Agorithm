from symmetries.objects.variables import Variable, Function, Power
from symmetries.objects.add import Add
from symmetries.objects.mul import Mul


x = Variable('x')
y = Variable('y')
p_of_x = Function('p', dependencies=(x,))
p_of_x_dx = Function('p', dependencies=(x,), derivatives=(x,))


def test_variable_x():
    assert x.__str__() == 'x'


def test_x_plus_x():
    assert (x+x).__str__() == '2x'


def test_x_plus_0():
    assert (x+0).__str__() == 'x'


def test_x_minus_x():
    assert (x-x).__str__() == '0'
    assert (x-x) == 0


def test_x_plus_2():
    assert (x+2).__str__() == 'x+2'
    assert isinstance(x+2, Add)


def test_x_times_2():
    assert (x*2).__str__() == '2x'
    assert isinstance(x*2, Mul)


def test_x_times_0():
    assert (x*0).__str__() == '0'
    assert x*0 == 0


def test_x_times_y():
    assert (x*y).__str__() == 'x*y'


def test_x_times_x():
    assert (x*x).__str__() == 'x^2'
    assert (x**2).__str__() == 'x^2'
    assert x**2 == x*x
    assert isinstance(x*x, Power)


def test_x_third_power():
    assert (x**3).__str__() == 'x^3'
    assert x**3 == x*x*x


def test_x_var_equals_2x():
    assert (x*2) == x


def test_x_var_not_equals_p():
    assert x != p_of_x


def test_function_p():
    assert p_of_x.__str__() == 'p(x)'


def test_p_times_x():
    assert (p_of_x*x).__str__() == 'p(x)*x'
    assert isinstance(p_of_x*x, Mul)


def test_x_times_x_second_power():
    assert ((x**2)**2).__str__() == 'x^4'


def test_x_powers():
    assert ((x**2)*(x**3)).__str__() == 'x^5'
    assert ((x**3)**2).__str__() == 'x^6'


def test_x_times_x_plus_1():
    assert (x*(x+1)).__str__() == 'x^2+x'


def test_x_minus_8():
    assert (x-8).__str__() == 'x-8'


def test_x_minus_8x():
    assert (x-8*x).__str__() == '-7x'


def test_x_add_8_plus_x():
    assert (x+(x+8)).__str__() == '2x+8'


def test_function_w_derivatives():
    assert p_of_x_dx.__str__() == 'p(x)_x'


def test_1_minus_x():
    assert (1-x).__str__() == '-x+1'


def test_1_plus_x():
    assert (1+x).__str__() == 'x+1'


def test_minus_add():
    assert (x-(y+y*x+2)).__str__() == '-y-x*y-2+x'

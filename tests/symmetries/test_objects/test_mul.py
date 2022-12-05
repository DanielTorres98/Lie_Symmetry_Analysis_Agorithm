from symmetries.objects.variables import Variable, Function, Power
from symmetries.objects.add import Add
from symmetries.objects.mul import Mul

x = Variable('x')
y = Variable('y')
p_of_x = Function('p', dependencies=(x,))

x_times_p = x*p_of_x


def test_x_times_minus_1():
    assert (x*-1).__str__() == '-x'
    assert isinstance(x*-1, Mul)


def test_int_times_mul():
    assert (3*x_times_p).__str__() == '3p(x)*x'


def test_x_times_2_squared():
    assert ((x*2)*(x*2)).__str__() == '4x^2'
    assert isinstance((x*2)*(x*2), Mul)


def test_x_times_p_squared():
    assert (x_times_p*x_times_p).__str__() == 'p(x)^2*x^2'
    assert isinstance(x_times_p*x_times_p, Mul)


def test_y_times_x_times_p():
    assert (y*x_times_p).__str__() == 'p(x)*x*y'
    assert isinstance(y*x_times_p, Mul)
    assert (y*x_times_p) == (x_times_p*y)


def test_y_times_x_times_x_plus_1():
    assert ((y*x)*(x+1)).__str__() == 'x^2*y+x*y'


def test_y_times_x_plus_8():
    assert ((y*x)+8).__str__() == 'x*y+8'
    assert ((x*y)+8).__str__() == 'x*y+8'


def test_y_times_x_plus_0():
    assert ((y*x)+0).__str__() == 'x*y'


def test_y_times_x_plus_x_plus_1():
    assert ((y*x)+(x+1)).__str__() == 'x+1+x*y'


def test_y_times_x_times_x():
    assert ((y*x)*x).__str__() == 'x^2*y'


def test_y_times_x_times_y():
    assert ((y*x)*y).__str__() == 'x*y^2'


def test_y_times_x_times_3():
    assert ((y*x)*3).__str__() == '3x*y'


def test_y_times_x_times_0():
    assert ((y*x)*0).__str__() == '0'
    assert ((y*x)*0) == 0


def test_y_times_x_then_times_p():
    assert ((y*x)*p_of_x).__str__() == 'p(x)*x*y'


def test_y_times_x_minus_same():
    assert ((y*x)+(-1*y*x)).__str__() == '0'
    assert ((y*x)+(-1*y*x)) == 0

def test_y_times_x_minus_x():
    assert ((y*x)+(-1*x)).__str__() == 'x*y+-x'
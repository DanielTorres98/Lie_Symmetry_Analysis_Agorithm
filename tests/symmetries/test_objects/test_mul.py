from symmetries.objects.variables import Variable, Function, Power
from symmetries.objects.add import Add
from symmetries.objects.mul import Mul

x = Variable('x')
y = Variable('y')
p_of_x = Function('p', dependencies=(x,))

x_times_p = x*p_of_x
x_times_2 = x*2
x_plus_p = x+p_of_x
x_plus_2 = x+2

def test_x_times_2_squared():
    assert (x_times_2*x_times_2).__str__()=='4x^2'
    assert isinstance(x_times_2*x_times_2, Mul)

def test_x_plus_2_squared():
    assert (x_plus_2*x_plus_2).__str__()=='x^2+4x+4'

def test_x_x_plus_p_squared():
    assert (x_plus_p*x_plus_p).__str__()=='x^2+2p(x)*x+p(x)^2'

def test_x_plus_2_times_x():
    assert (x_plus_2*x).__str__()=='x^2+2x'
    assert isinstance(x_plus_2*x, Add)

def test_x_times_p_squared():
    assert (x_times_p*x_times_p).__str__()=='p(x)^2*x^2'
    assert isinstance(x_times_p*x_times_p, Mul)

def test_y_times_x_times_p():
    assert (y*x_times_p).__str__()=='p(x)*x*y'
    assert isinstance(y*x_times_p, Mul)
    assert (y*x_times_p)==(x_times_p*y)

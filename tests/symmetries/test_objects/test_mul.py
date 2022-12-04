from symmetries.objects.variables import Variable, Function, Power
from symmetries.objects.add import Add
from symmetries.objects.mul import Mul

x = Variable('x')
y = Variable('y')
p_of_x = Function('p', dependencies=(x,))

x_plus_y = x+y
x_squared = x**2
x_times_2 = x*2
x_plus_2 = x+2

def test_x_times_2_squared():
    assert x_times_2*x_times_2==x_times_2**2
    assert x_times_2*x_times_2.__str__()=='4x^2'
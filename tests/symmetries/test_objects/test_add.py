from symmetries.objects.variables import Variable, Function, Power
from symmetries.objects.add import Add

x = Variable('x')
y = Variable('y')
p_of_x = Function('p', dependencies=(x,))


x_plus_p = x+p_of_x
x_plus_2 = x+2
x_plus_y = x+y


def test_x_plus_2_squared():
    assert (x_plus_2*x_plus_2).__str__() == 'x^2+4x+4'


def test_x_x_plus_p_squared():
    assert (x_plus_p*x_plus_p).__str__() == 'x^2+2p(x)*x+p(x)^2'


def test_x_plus_2_times_x():
    assert (x_plus_2*x).__str__() == 'x^2+2x'
    assert isinstance(x_plus_2*x, Add)


def test_x_plus_2_times_2x():
    assert ((x+2)*(2*x)).__str__() == '2x^2+4x'
    assert isinstance((x+2)*(2*x), Add)


def test_x_2_plus_x_2():
    assert ((x+2)+(x+2)).__str__() == '2x+4'


def test_x_2_times_4():
    assert ((x+2)*4).__str__() == '4x+8'
    assert (4*(x+2)).__str__() == '4x+8'


def test_x_2_times_0():
    assert ((x+2)*0).__str__() == '0'
    assert ((x+2)*0) == 0


def test_1_minus_add():
    assert (1-(y+y*x+2)).__str__() == '-y-x*y-1'


def test_1_plus_add():
    assert (1+(y+y*x+2)).__str__() == 'y+x*y+3'


def test_add_minus_add():
    assert ((x*9-y*7)-(y+y*x+2)).__str__() == '9x-8y-x*y-2'


def test_add_plus_add():
    assert ((8*x+y-9)+(y-x+4)).__str__() == '7x+2y-5'


def test_add_remove_var_equals_0():
    assert ((8*x+y-9)-(y-x+4)).__str__() == '9x-13'


def test_add_power_5():
    assert ((1-x)**5).__str__() == '-x^5+5x^4-10x^3+10x^2-5x+1'
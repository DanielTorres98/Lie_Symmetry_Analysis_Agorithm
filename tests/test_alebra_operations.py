from utils import algebra

def test_is_zero(is_zero_data_test):
    zero_term = {"coefficient": 4,
                 "constants": [2, 3],
                 "derivatives": [0, 1, 3, 2],
                 "variable": 'xse'}
    term = is_zero_data_test
    assert term["answer"] == algebra.is_zero(zero_term, term)

def test_drop_constants(drop_constants_data_test):
    eqn = drop_constants_data_test
    assert eqn["result"] == algebra.drop_constants(eqn["test"])
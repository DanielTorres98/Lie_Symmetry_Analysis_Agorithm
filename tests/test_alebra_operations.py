from Notebooks.utils import algebra

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

def test_compare_derivatives(compare_derivatives_data_test):
    D1 = compare_derivatives_data_test["zero_term"]
    D2 = compare_derivatives_data_test["term"]
    assert compare_derivatives_data_test["result"] == algebra.compare_derivatives(D1, D2)

def test_key_ordering(key_ordering_data_test):
    original_list = key_ordering_data_test["original"]
    result = key_ordering_data_test["result"]
    ordered_list = algebra.key_ordering(original_list)
    assert result == ordered_list
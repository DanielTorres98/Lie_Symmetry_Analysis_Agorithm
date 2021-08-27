import json

from pytest import fixture


data_is_zero = "data_test/data_algebra/data_is_zero.json"

data_drop_constants = "data_test/data_algebra/data_drop_constants.json"

data_compare_derivatives = "data_test/data_algebra/data_compare_derivatives.json"

def load_test_data(path):
    with open(path) as data_file:
        data = json.load(data_file)
        return data

@fixture(params=load_test_data(data_is_zero))
def is_zero_data_test(request):
    data = request.param
    return data

@fixture(params=load_test_data(data_drop_constants))
def drop_constants_data_test(request):
    data = request.param
    return data

@fixture(params=load_test_data(data_compare_derivatives))
def compare_derivatives_data_test(request):
    data = request.param
    return data
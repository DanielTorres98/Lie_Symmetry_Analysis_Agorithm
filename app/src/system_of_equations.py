"""This module contains the AbstractBaseClass for the systems of equations."""
from abc import ABC


class SystemOfEquations(ABC):

    def __init__(self, independent_variables: list, dependent_variables: list, constants: list):
        self.independent_variables = independent_variables
        self.dependent_variables = dependent_variables
        self.all_variables = independent_variables + dependent_variables
        self.constants = constants

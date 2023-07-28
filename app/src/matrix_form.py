
"""Module containing functions for turning the equations into matrix form"""

class MatrixForm():
    # latex_form = latex.format_determining_equations(False)
    def __init__(self, latex_form) -> None:
        self.equations = self._split(latex_form)

    def _split(self, latex_form):
        equations = latex_form.split("\\\n")
        equations = [eq.split('=')[0] for eq in equations]
        return [eq.replace('\\', '').split('+') for eq in equations]

    def get_variables():
        pass
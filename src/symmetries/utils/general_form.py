
from copy import deepcopy
import regex as re
from symmetries.objects.system import System
from symmetries.objects.determining_equations import DeterminingEquations
from symmetries.utils.symbolic import sym_det_eqn

class GeneralForm():
    def __init__(self, system: System, model:DeterminingEquations)  -> None:
        self.system = system
        self.model = model
        self.general_form = dict

    def obtain_general_form(self):
        general_form = {}
        infinitesimals = self.model.infinitesimals_ind + self.model.infinitesimals_dep

        for inft in infinitesimals:
            l = re.split(r'\W+', str(inft), len(self.model.all_variables)+1)
            if l[0] not in general_form:
                general_form[l[0]] = {l[1]: deepcopy(self.model.all_variables)}
            else:
                general_form[l[0]][l[1]] = deepcopy(self.model.all_variables)

        self.general_form = general_form

    def find_first_derivative_equals_0(self):
        to_delete = []
        for k, v in self.system.determining_equations.items():
            if len(v)==1:
                l = v[0]
                if sum(l['derivatives'])==1:
                    var = [self.model.all_variables[n] for n, d in enumerate(l['derivatives']) if d][0]
                    self.general_form[l['variable'][:-1]][l['variable'][-1]].remove(var)
                    to_delete.append(k)

        for k in to_delete:
            del self.system.determining_equations[k]

    def print_matrix(self):
        return sym_det_eqn(
            self.system.determining_equations, 
            self.model.independent_variables, 
            self.model.dependent_variables, 
            self.model.constants)
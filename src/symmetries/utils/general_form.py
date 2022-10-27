
from copy import deepcopy
import regex as re
from symmetries.objects.system import System
from symmetries.objects.determining_equations import DeterminingEquations
from symmetries.utils.symbolic import sym_det_eqn


class GeneralForm():
    def __init__(self, system: System, model: DeterminingEquations) -> None:
        self.model = model
        self.general_form = self.obtain_general_form()
        self.determining_equations = deepcopy(system.determining_equations)
        self.deleted = {}

    def obtain_general_form(self):
        general_form = {}
        infinitesimals = self.model.infinitesimals_ind + self.model.infinitesimals_dep

        for inft in infinitesimals:
            l = re.split(r'\W+', str(inft), len(self.model.all_variables)+1)
            if l[0] not in general_form:
                general_form[l[0]] = {l[1]: deepcopy(self.model.all_variables)}
            else:
                general_form[l[0]][l[1]] = deepcopy(self.model.all_variables)

        return general_form

    def find_first_derivative_equals_0(self):
        for k, v in self.determining_equations.items():
            if len(v) == 1:
                l = v[0]
                if sum(l['derivatives']) == 1:
                    var = [self.model.all_variables[n]
                        for n, d in enumerate(l['derivatives']) if d][0]
                    if not (l['variable'] in self.deleted and var in self.deleted[l['variable']]):
                        self.general_form[l['variable'][:-1]
                                        ][l['variable'][-1]].remove(var)
                        if l['variable'] not in self.deleted:
                            self.deleted[l['variable']] = [var]
                            print('deleting', l['variable'], var)
                        else:
                            self.deleted[l['variable']].append(var)


    def find_deleted_items_in_equations(self):
        for k, v in self.determining_equations.items():
            if len(v) > 1:
                values = deepcopy(v)
                for item in values:
                    if item['variable'] in self.deleted:
                        var = [self.model.all_variables[n]
                                for n, d in enumerate(item['derivatives']) if d][0]
                        if var in self.deleted[item['variable']]:
                            print('found', var, 'in', item['variable'], 'eq', k)
                            v.remove(item)

    def delete_second_derivatives():
        pass
                            
    def print_matrix(self):
        print('general form:', self.general_form)
        print('already deleted:', self.deleted)
        return sym_det_eqn(
            self.determining_equations,
            self.model.independent_variables,
            self.model.dependent_variables,
            self.model.constants)

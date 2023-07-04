# BilinearProblem Setup
# Xiaoyi Gu


import numpy as np
from math import sqrt  # for readability
from itertools import chain, combinations

import warnings
from typing import Union, Optional

import gurobipy as gp
from gurobipy import GRB
import mosek.fusion as MF
# import pyomo.environ as pyo


EPS_SAFE = 1e-4


def powerset(iterable):
    # courtesy https://docs.python.org/3/library/itertools.html
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))


class BilinearGurobiModel:
    def __init__(self, name=None, OutputFlag=1):
        self.OutputFlag = OutputFlag
        self.m = None  # model
        self.x = None
        self.y = None
        self.w = None  # xy = w
        self.v = None  # xy = v^2
        # w <= x, w <= y, w >= x+y-1, w >= 0
        # v <= sqrt(xy)
        # ADD McCormick w = v^2
        # w <= min(x,y) <= v IMPLIED
        # w >= x+y-1 >= 2v-1 ?
        # MIGHT ADD w >= v^2
        self.u = None  # min{x,y} = u
        self.b_lhs = None  # bilinear lhs
        self.obj = None

        self._init_grb_model(name)

    def _init_grb_model(self, name: str = ""):
        if type(name) is str:
            self.m = gp.Model(name)
        else:
            self.m = gp.Model("")
        self.m.setParam("OutputFlag", self.OutputFlag)

    def reset(self, name: str = ""):
        self.__init__(name)


class BilinearMinimalCover:
    def __init__(self, bi_index: int = 0):
        self.label_i = -1
        self.label_a = 2
        self.label_0 = 0
        self.label_1 = 1
        self.label = None  # label 2, 0, 1, -1 for active, 0, 1, inactive (not present)
        self.a_i_0 = None  # smallest a_i_0 > delta if exists
        self.delta = None  # sum a_i_active - d
        self.rhs = None  # d
        self.lp = None  # l_+
        self.ln = None  # l_-

        self.grb_var = None  # grb var
        self.grb_lhs = None  # grb lhs (i.e. sum of var)

        self.bi_index = bi_index


class BilinearProblem:
    def __init__(self, name=None, OutputFlag=1):
        self.name = name
        self.size = 0

        self.x_cost = []
        self.y_cost = []
        self.bi_costs = []
        self.bi_rhss = []
        self.bi_size = 0

        self.min_cover_list = []

        self.grb = None

        self.OutputFlag = OutputFlag

        self.test_grb = None

    def grb_init(self):
        self.grb = BilinearGurobiModel(self.name, OutputFlag=self.OutputFlag)

    def add_costs(self, x_cost: list, y_cost: list):
        if len(x_cost) != len(y_cost):
            raise RuntimeError("x_cost, y_cost - Length NOT Match!")
        self.size = len(x_cost)
        self.x_cost = x_cost
        self.y_cost = y_cost

    def add_bilinear(self, bi_cost: list, bi_rhs: float = 0.0):
        if len(bi_cost) != self.size:
            raise RuntimeError("bi_cost - Length NOT Match!")
        self.bi_costs.append(bi_cost)
        self.bi_rhss.append(bi_rhs)
        self.bi_size += 1

    def add_bilinears(self, bi_costs: list, bi_rhss: list = 0.0):
        if type(bi_rhss) is list:
            if len(bi_costs) != len(bi_rhss):
                raise RuntimeError("bi_cost, bi_rhs - Lengths NOT Match!")
            for i in range(len(bi_rhss)):
                self.add_bilinear(bi_costs[i], bi_rhss[i])
        else:
            for i in range(len(bi_rhss)):
                self.add_bilinear(bi_costs[i], bi_rhss)

    def grb_reset(self):
        self.grb.reset(self.name)

    def grb_init_basic(self):
        self.grb_init()
        self._grb_set_var()
        self._grb_set_obj()
        self._grb_set_McCormick()
        self._grb_set_bilinear()
        self._grb_set_min()
        self._grb_set_SOS()

    def add_all_minimal_cover_cuts_grb(self, debug: bool = False):
        # add all valid minimal cover cuts with only one active element within minimal cover
        iter_powerset = list(powerset(range(self.size)))
        for itr in range(1, 2 ** self.size):
            for bi_index in range(self.bi_size):
                cover_element_index = list(iter_powerset[itr])
                tmp_bi_cost = self.bi_costs[bi_index]
                tmp_bi_rhs = self.bi_rhss[bi_index]
                fixed_sum = sum(tmp_bi_cost[ind] for ind in cover_element_index)
                delta = fixed_sum - tmp_bi_rhs
                if delta < 1e-4:
                    continue
                for ind_active in cover_element_index:
                    if tmp_bi_cost[ind_active] < delta:
                        continue
                    label = [(ind in cover_element_index) + (ind == ind_active) for ind in range(self.size)]
                    if self.check_minimal_cover(label, bi_index) == 0:
                        self.add_minimal_cover_grb_checked(label, bi_index)

    def add_all_minimal_cover_cuts_range_grb(self, iterator, debug: bool = False):
        # add all valid minimal cover cuts for an iterator with only one active element within minimal cover
        tmp_size = len(iterator)
        iter_powerset = list(powerset(iterator))
        for itr in range(1, 2 ** tmp_size):
            for bi_index in range(self.bi_size):
                cover_element_index = list(iter_powerset[itr])
                tmp_bi_cost = self.bi_costs[bi_index]
                tmp_bi_rhs = self.bi_rhss[bi_index]
                fixed_sum = sum(tmp_bi_cost[ind] for ind in cover_element_index)
                delta = fixed_sum - tmp_bi_rhs
                if delta < 1e-4:
                    continue
                for ind_active in cover_element_index:
                    if tmp_bi_cost[ind_active] < delta:
                        continue
                    label = [(ind in cover_element_index) + (ind == ind_active) for ind in range(self.size)]
                    if self.check_minimal_cover(label, bi_index) == 0:
                        self.add_minimal_cover_grb_checked(label, bi_index)

    def add_minimal_cover_cuts_enumeration_grb(self, bi_index, iterator, debug: bool = False):
        # add all valid minimal cover cuts for an iterator with only one active element within minimal cover
        tmp_size = len(iterator)
        iter_powerset = list(powerset(iterator))
        if bi_index >= self.bi_size:
            raise RuntimeError("add_minimal_cover_cuts_enumeration_grb ERROR! bi_index out of range!")

        for itr in range(1, 2 ** tmp_size):
            cover_element_index = list(iter_powerset[itr])
            tmp_bi_cost = self.bi_costs[bi_index]
            tmp_bi_rhs = self.bi_rhss[bi_index]
            fixed_sum = sum(tmp_bi_cost[ind] for ind in cover_element_index)
            delta = fixed_sum - tmp_bi_rhs
            if delta < 1e-4:
                continue
            for ind_active in cover_element_index:
                if tmp_bi_cost[ind_active] < delta:
                    continue
                label = [(ind in cover_element_index) + (ind == ind_active) for ind in range(self.size)]
                if debug:
                    self.add_minimal_cover_grb(label, bi_index)
                else:
                    self.add_minimal_cover_grb_checked(label, bi_index)

    def add_minimal_cover_grb(self, label: list, bi_index: int = 0):
        # add minimal cover with its cut
        # input - label: label for minimal cover
        #         bi_index: index for bilinear constraint
        mc = BilinearMinimalCover(bi_index=bi_index)
        tmp_bi_cost = self.bi_costs[bi_index]
        tmp_bi_rhs = self.bi_rhss[bi_index]
        mc.label = label
        fixed_sum = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_1)
        mc.rhs = tmp_bi_rhs - fixed_sum
        if mc.rhs < 1e-4:
            error_string = "Minimal Cover RHS Too Small! rhs =" + str(mc.rhs) + "< 1e-4"
            raise RuntimeError(error_string)
        mc.delta = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_a) - mc.rhs
        if mc.delta < 1e-4:
            error_string = "Minimal Cover Delta Too Small! delta =" + str(mc.delta) + "< 1e-4"
            raise RuntimeError(error_string)
        if [ind for ind, val in enumerate(mc.label) if val == mc.label_a and tmp_bi_cost[ind] < mc.delta - 1e-4]:
            raise RuntimeError("Minimal Cover Cost Smaller Than Delta!")
        mc.ln = 1 / mc.delta
        a_i_0_candidate = list(filter(lambda x: x >= mc.delta + 1e-4, tmp_bi_cost))
        if a_i_0_candidate:
            mc.a_i_0 = min(a_i_0_candidate)
            mc.lp = (sqrt(mc.a_i_0 / (mc.a_i_0 - mc.delta)) + 1) / mc.delta
        else:
            mc.lp = 1 / mc.delta
        self.min_cover_list.append(mc)
        self._grb_add_valid(mc)

    def check_minimal_cover(self, label: list, bi_index: int = 0):
        # check if the proposed set forms a minimal cover
        # input - label: label for minimal cover
        #         bi_index: index for bilinear constraint
        # return -  -1: rhs too small
        #           -2: minimal cover delta too small
        #           1: minimal cover cost smaller than delta
        #           0: no problem

        return self.check_minimal_cover_generate(label, bi_index)[0]

    def check_minimal_cover_generate(self, label: list, bi_index: int = 0):
        # check if the proposed set forms a minimal cover
        # input - label: label for minimal cover
        #         bi_index: index for bilinear constraint
        # return -  -1: rhs too small
        #           -2: minimal cover delta too small
        #           1: minimal cover cost smaller than delta
        #           0: no problem
        #        -  mc: generated mc
        #        -  None: not a minimal cover
        mc = BilinearMinimalCover(bi_index=bi_index)
        tmp_bi_cost = self.bi_costs[bi_index]
        tmp_bi_rhs = self.bi_rhss[bi_index]
        mc.label = label
        fixed_sum = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_1)
        mc.rhs = tmp_bi_rhs - fixed_sum
        if mc.rhs < 1e-4:
            return -1, None  # rhs too small

        mc.delta = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_a) - mc.rhs

        if mc.delta < 1e-4:
            return -2, None  # minimal cover delta too small

        if [ind for ind, val in enumerate(mc.label) if val == mc.label_a and tmp_bi_cost[ind] < mc.delta - 1e-4]:
            return 1, None  # minimal cover cost smaller than delta

        mc.ln = 1 / mc.delta
        a_i_0_candidate = list(filter(lambda x: x >= mc.delta + 1e-4, tmp_bi_cost))
        if a_i_0_candidate:
            mc.a_i_0 = min(a_i_0_candidate)
            mc.lp = (sqrt(mc.a_i_0 / (mc.a_i_0 - mc.delta)) + 1) / mc.delta
        else:
            mc.lp = 1 / mc.delta

        return 0, mc

    def add_minimal_cover_grb_generated(self, mc: Optional[BilinearMinimalCover]):
        # add a checked minimal cover with its cut
        # input - mc: BilinearMinimalCover
        if mc is None:
            warnings.warn("None input as BilinearMInimalCover!")
            return

        self.min_cover_list.append(mc)
        self._grb_add_valid(mc)

    def add_minimal_cover_grb_checked(self, label: list, bi_index: int = 0):
        # add a checked minimal cover with its cut
        # omit the steps of checking validity
        # input - label: label for minimal cover
        mc = BilinearMinimalCover(bi_index=bi_index)
        tmp_bi_cost = self.bi_costs[bi_index]
        tmp_bi_rhs = self.bi_rhss[bi_index]

        mc.label = label
        fixed_sum = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_1)
        mc.rhs = tmp_bi_rhs - fixed_sum
        mc.delta = sum(tmp_bi_cost[ind] for ind, val in enumerate(mc.label) if val == mc.label_a) - mc.rhs
        mc.ln = 1 / mc.delta
        a_i_0_candidate = list(filter(lambda x: x >= mc.delta + 1e-4, tmp_bi_cost))
        if a_i_0_candidate:
            mc.a_i_0 = min(a_i_0_candidate)
            mc.lp = (sqrt(mc.a_i_0 / (mc.a_i_0 - mc.delta)) + 1) / mc.delta
        else:
            mc.lp = 1 / mc.delta

        # self.min_cover_list.append(mc)
        # self._grb_add_valid(mc)
        self.add_minimal_cover_grb_generated(mc)

    def _grb_set_var(self):
        # initiate gurobi variables
        # x, y, w (McCormick xy), v (v^2 <= xy), u (u <= min{x,y})
        self.grb.x = self.grb.m.addVars(self.size, lb=0, ub=1, name="x")
        self.grb.y = self.grb.m.addVars(self.size, lb=0, ub=1, name="y")
        self.grb.w = self.grb.m.addVars(self.size, lb=0, ub=1, name="w")
        self.grb.v = self.grb.m.addVars(self.size, lb=0, ub=1, name="v")
        self.grb.u = self.grb.m.addVars(self.size, lb=0, ub=1, name="u")
        self.grb.m.update()

    def _grb_set_obj(self):
        # initiate gurobi objective
        self.grb.obj = sum(self.grb.x[j] * self.x_cost[j] for j in range(self.size)) + \
                       sum(self.grb.y[j] * self.y_cost[j] for j in range(self.size))
        self.grb.m.setObjective(self.grb.obj, GRB.MINIMIZE)
        self.grb.m.update()

    def _grb_set_SOS(self):
        # initiate v: gurobi SOS: v^2 <= xy
        self.grb.m.addConstrs(self.grb.v[j] * self.grb.v[j] <= self.grb.x[j] * self.grb.y[j]
                              for j in range(self.size))
        self.grb.m.update()

    def _grb_set_min(self):
        # initiate u: gurobi min: u <= x,y
        self.grb.m.addConstrs(self.grb.x[j] >= self.grb.u[j] for j in range(self.size))
        self.grb.m.addConstrs(self.grb.y[j] >= self.grb.u[j] for j in range(self.size))
        self.grb.m.update()

    def _grb_set_McCormick(self):
        # initiate w: McCormick xy
        # xy >= 0 -- unnecessary
        # xy >= x+y-1
        self.grb.m.addConstrs(self.grb.w[j] >= self.grb.x[j] + self.grb.y[j] - 1
                              for j in range(self.size))
        # xy <= x
        self.grb.m.addConstrs(self.grb.w[j] <= self.grb.x[j] for j in range(self.size))
        # xy <= y
        self.grb.m.addConstrs(self.grb.w[j] <= self.grb.y[j] for j in range(self.size))
        self.grb.m.update()

    def _grb_set_bilinear(self):
        # initiate bilinear with w (McCormick xy)
        self.grb.b_lhs = [None for i in range(len(self.bi_costs))]
        for bi_index in range(self.bi_size):
            tmp_bi_cost = self.bi_costs[bi_index]
            tmp_bi_rhs = self.bi_rhss[bi_index]
            self.grb.b_lhs[bi_index] = sum(self.grb.w[j] * tmp_bi_cost[j] for j in range(self.size))
            self.grb.m.addConstr(self.grb.b_lhs[bi_index] >= tmp_bi_rhs)
        self.grb.m.update()

    def _grb_check_valid(self, mc: BilinearMinimalCover):
        # check if the bilinear minimal cover cut is violated by the current solution
        # return the exact value (negative means violated)
        tmp_var = [0.0] * self.size
        tmp_bi_cost = self.bi_costs[mc.bi_index]
        tmp_bi_rhs = self.bi_rhss[mc.bi_index]
        tmp_x = self.grb_get_x()
        tmp_y = self.grb_get_y()

        for j in range(self.size):
            if mc.label[j] == mc.label_i:  # inactive (not present)
                if tmp_bi_cost[j] != 0:
                    raise RuntimeError(f"Minimal Cover Label Incorrectly Deactivated!")
                continue
            elif mc.label[j] == mc.label_a:  # active
                tmp = tmp_bi_cost[j] - mc.delta  # d_j
                if tmp < 0:
                    tmp = 0
                tmp_var[j] = (tmp_bi_cost[j] + sqrt(tmp_bi_cost[j] * tmp)) / mc.delta * (
                        sqrt(tmp_x[j] * tmp_y[j]) - 1)
            elif mc.label[j] == mc.label_0:  # 0
                if tmp_bi_cost[j] >= 0:  # 0 pos
                    tmp_var[j] = mc.lp * tmp_bi_cost[j] * min(tmp_x[j], tmp_y[j])
                else:  # 0 neg
                    tmp_var[j] = min(
                        - mc.ln * tmp_bi_cost[j] * (1 - tmp_x[j] - tmp_y[j]),
                        - mc.lp * tmp_bi_cost[j] * (1 - tmp_x[j] - tmp_y[j]) + mc.lp * mc.delta - 1,
                        0
                    )
            elif mc.label[j] == mc.label_1:  # 1
                if tmp_bi_cost[j] >= 0:  # 1 pos
                    if mc.a_i_0 is None or tmp_bi_cost[j] < mc.a_i_0:
                        tmp_var[j] = min(
                            mc.lp * tmp_bi_cost[j] * (min(tmp_x[j], tmp_y[j]) - 1) + mc.lp * mc.delta - 1,
                            mc.ln * tmp_bi_cost[j] * (min(tmp_x[j], tmp_y[j]) - 1)
                        )
                    else:
                        tmp = tmp_bi_cost[j] - mc.delta  # d_j
                        tmp_var[j] = min(
                            mc.lp * tmp_bi_cost[j] * (min(tmp_x[j], tmp_y[j]) - 1) + mc.lp * mc.delta - 1,
                            mc.ln * tmp_bi_cost[j] * (min(tmp_x[j], tmp_y[j]) - 1),
                            mc.lp * sqrt(tmp * tmp_bi_cost[j] * tmp_x[j] * tmp_y[j]) - mc.lp * tmp - 1,
                            (tmp_bi_cost[j] + sqrt(tmp_bi_cost[j] * tmp)) / mc.delta * (sqrt(tmp_x[j] * tmp_y[j]) - 1)
                        )
                else:  # 1 neg
                    tmp_var[j] = - mc.lp * tmp_bi_cost[j] * min(2 - tmp_x[j] - tmp_y[j], 1)
            else:
                raise RuntimeError(f"Unknown Minimal Cover Label! Label {mc.label[j]} is illegal.")
        return sum(tmp_var) + 1 + EPS_SAFE

    def _grb_add_valid(self, mc: BilinearMinimalCover):
        # add valid cut to gurobi model given a minimal cover
        tmp_index_list = [j for j in range(self.size) if mc.label[j] != mc.label_i]   # remove inactive index
        tmp_len = len(tmp_index_list)
        mc.grb_var = self.grb.m.addVars(tmp_len, lb=-float('inf'), ub=float('inf'))

        tmp_bi_cost = self.bi_costs[mc.bi_index]
        tmp_bi_rhs = self.bi_rhss[mc.bi_index]

        for ind in range(tmp_len):
            j = tmp_index_list[ind]
            if mc.label[j] == mc.label_a:  # active
                tmp = tmp_bi_cost[j] - mc.delta  # d_j
                if tmp < 0:
                    tmp = 0
                self.grb.m.addConstr(mc.grb_var[ind] ==
                                     (tmp_bi_cost[j] + sqrt(tmp_bi_cost[j] * tmp)) / mc.delta * (self.grb.v[j] - 1)
                                     )
            elif mc.label[j] == mc.label_0:  # 0
                if tmp_bi_cost[j] >= 0:  # 0 pos
                    self.grb.m.addConstr(mc.grb_var[ind] == mc.lp * tmp_bi_cost[j] * self.grb.u[j])
                else:  # 0 neg
                    self.grb.m.addConstr(
                        mc.grb_var[ind] <= - mc.ln * tmp_bi_cost[j] * (1 - self.grb.x[j] - self.grb.y[j]))
                    self.grb.m.addConstr(
                        mc.grb_var[ind] <= - mc.lp * tmp_bi_cost[j] * (1 - self.grb.x[j] - self.grb.y[j])
                        + mc.lp * mc.delta - 1)
                    mc.grb_var[ind].ub = 0
            elif mc.label[j] == mc.label_1:  # 1
                if tmp_bi_cost[j] >= 0:  # 1 pos
                    if mc.a_i_0 is None or tmp_bi_cost[j] < mc.a_i_0:
                        self.grb.m.addConstr(mc.grb_var[ind] <= mc.lp * tmp_bi_cost[j] * (self.grb.u[j] - 1)
                                             + mc.lp * mc.delta - 1)
                        self.grb.m.addConstr(mc.grb_var[ind] <= mc.ln * tmp_bi_cost[j] * (self.grb.u[j] - 1))
                    else:
                        tmp = tmp_bi_cost[j] - mc.delta  # d_j
                        self.grb.m.addConstr(mc.grb_var[ind] <= mc.lp * tmp_bi_cost[j] * (self.grb.u[j] - 1)
                                             + mc.lp * mc.delta - 1)
                        self.grb.m.addConstr(mc.grb_var[ind] <= mc.ln * tmp_bi_cost[j] * self.grb.u[j] - 1)
                        self.grb.m.addConstr(mc.grb_var[ind] <= sqrt(tmp_bi_cost[j] * tmp) * mc.lp * self.grb.v[j]
                                             - mc.lp * tmp - 1)
                        self.grb.m.addConstr(mc.grb_var[ind] <=
                                             (tmp_bi_cost[j] + sqrt(tmp_bi_cost[j] * tmp))
                                             / mc.delta * (self.grb.v[j] - 1)
                                             )
                else:  # 1 neg
                    mc.grb_var[ind].ub = - mc.lp * tmp_bi_cost[j]
                    self.grb.m.addConstr(
                        mc.grb_var[ind] <= - mc.lp * tmp_bi_cost[j] * (2 - self.grb.x[j] - self.grb.y[j]))
            else:
                raise RuntimeError(f"Unknown Minimal Cover Label! Label {mc.label[j]} is illegal.")
        mc.grb_lhs = sum(mc.grb_var[ind] for ind in range(tmp_len)) + 1 + EPS_SAFE
        self.grb.m.addConstr(mc.grb_lhs >= 0)
        self.grb.m.update()

    def grb_opt(self):
        self.grb.m.update()
        self.grb.m.optimize()

    def _test_grb_bilinear_opt_tmp(self, OutputFlag=1, extra_print=0):
        # solve the nonconvex bilinear problem using gurobi builtin nonconvex bilinear solver
        # destroy the instance!
        self.grb.m.params.NonConvex = 2
        self.grb.m.addConstrs(self.grb.w[j] == self.grb.x[j] * self.grb.y[j] for j in range(self.size))
        self.grb.m.update()

        self.grb.m.optimize()

    # def test_pyomo_bilinear_opt(self, solver="baron", extra_print=0):
    #     # solve the nonconvex bilinear problem using BARON solver
    #     model = pyo.ConcreteModel()
    #     model.x = pyo.Var(range(self.size), within=pyo.NonNegativeReals, bounds=(0, 1), initialize=0)
    #     model.y = pyo.Var(range(self.size), within=pyo.NonNegativeReals, bounds=(0, 1), initialize=0)
    #     model.w = pyo.Var(range(self.size), within=pyo.NonNegativeReals, bounds=(0, 1), initialize=0)

    #     model.obj = pyo.Objective(expr=(
    #         sum(model.x[j] * self.x_cost[j] + model.y[j] * self.y_cost[j] for j in range(self.size))
    #     ), sense=pyo.minimize)

    #     model.cw = pyo.ConstraintList()
    #     for j in range(self.size):
    #         model.cw.add(model.w[j] == model.x[j] * model.y[j])

    #     model.cb = pyo.ConstraintList()
    #     for i in range(self.bi_size):
    #         if all(bi_cost == 0 for bi_cost in self.bi_costs[i]):
    #             continue
    #         model.cb.add(sum(model.w[j] * self.bi_costs[i][j] for j in range(self.size))
    #                      >= self.bi_rhss[i])

    #     model.display()
    #     if solver == None:
    #         return model
    #     elif solver == "gurobi":
    #         solver = pyo.SolverFactory(solver, solver_io="python")
    #         solver.options["NonConvex"] = 2
    #     else:
    #         solver = pyo.SolverFactory(solver)
    #     results = solver.solve(model)

    #     x_sol = [model.x[j].value for j in range(self.size)]
    #     y_sol = [model.y[j].value for j in range(self.size)]
    #     w_sol = [model.w[j].value for j in range(self.size)]

    #     real_lhs = [sum(x_sol[j] * y_sol[j] * self.bi_costs[i][j] for j in range(self.size))
    #                 for i in range(self.bi_size)]

    #     gap = [self.bi_rhss[i] - real_lhs[i] for i in range(self.bi_size)]

    #     if extra_print == 1:
    #         with np.printoptions(precision=3, suppress=True):
    #             print("gurobi nonconvex bilinear:")
    #             print("obj:")
    #             print(np.array(model.obj()))
    #             print("x:")
    #             print(np.array(x_sol))
    #             print("y:")
    #             print(np.array(y_sol))
    #             print("w:")
    #             print(np.array(w_sol))
    #             print("gap:")
    #             print(np.array(gap))

    #     return model.obj(), model, results

    def test_grb_opt_current(self, OutputFlag=1, extra_print=0, time_limit=120):
        self.grb.m.params.NonConvex = 2
        self.grb.m.params.TimeLimit = time_limit
        self.grb.m.addConstrs(
            self.grb.w[j] == self.grb.x[j] * self.grb.y[j]
            for j in range(self.size)
        )
        self.grb.m.update()

        self.grb.m.optimize()

    def test_grb_bilinear_opt(self, OutputFlag=1, extra_print=0, time_limit=120):
        # solve the nonconvex bilinear problem using gurobi builtin nonconvex bilinear solver
        model = gp.Model("")
        model.params.TimeLimit = time_limit

        model.setParam("OutputFlag", OutputFlag)
        model.params.NonConvex = 2
        x = model.addVars(self.size, lb=0, ub=1, name="x")
        y = model.addVars(self.size, lb=0, ub=1, name="x")
        w = model.addVars(self.size, lb=0, ub=1, name="w")
        model.update()

        obj = sum(x[j] * self.x_cost[j] for j in range(self.size)) + \
              sum(y[j] * self.y_cost[j] for j in range(self.size))
        model.setObjective(obj, GRB.MINIMIZE)

        model.addConstrs(w[j] == x[j] * y[j] for j in range(self.size))
        model.update()

        lhs = [sum(w[j] * self.bi_costs[i][j] for j in range(self.size)) for i in range(self.bi_size)]
        model.addConstrs(lhs[i] >= self.bi_rhss[i] for i in range(self.bi_size))
        model.update()

        model.optimize()

        self.test_grb = model

        x_sol = [x[j].X for j in range(self.size)]

        y_sol = [y[j].X for j in range(self.size)]

        w_sol = [w[j].X for j in range(self.size)]

        real_lhs = [sum(x_sol[j] * y_sol[j] * self.bi_costs[i][j] for j in range(self.size))
                    for i in range(self.bi_size)]

        gap = [self.bi_rhss[i] - real_lhs[i] for i in range(self.bi_size)]

        if extra_print == 1:
            with np.printoptions(precision=3, suppress=True):
                print("gurobi nonconvex bilinear:")
                print("obj:")
                print(np.array(model.objVal))
                print("x:")
                print(np.array(x_sol))
                print("y:")
                print(np.array(y_sol))
                print("w:")
                print(np.array(w_sol))
                print("bi_lhs:")
                print(np.array([_.getValue() for _ in lhs]))
                print("gap:")
                print(np.array(gap))

        out_objVal, out_objBound = model.objVal, model.objBound
        del model
        return out_objVal, out_objBound

    def test_grb_bilinear_opt_Richard(self, OutputFlag=1, extra_print=0, time_limit=120):
        # tighten root relaxation with Richard results
        # solve the nonconvex bilinear problem using gurobi builtin nonconvex bilinear solver
        model = gp.Model("")
        model.params.TimeLimit = time_limit

        model.setParam("OutputFlag", OutputFlag)
        x = model.addVars(self.size, lb=0, ub=1, name="x")
        y = model.addVars(self.size, lb=0, ub=1, name="x")
        w = model.addVars(self.size, lb=0, ub=1, name="w")
        v = model.addVars(self.size, lb=0, ub=1, name="w")
        model.update()

        obj = sum(x[j] * self.x_cost[j] for j in range(self.size)) + \
              sum(y[j] * self.y_cost[j] for j in range(self.size))
        model.setObjective(obj, GRB.MINIMIZE)

        # McCormick
        # xy >= x+y-1
        model.addConstrs(w[j] >= x[j] + y[j] - 1 for j in range(self.size))
        # xy <= x
        model.addConstrs(w[j] <= x[j] for j in range(self.size))
        # xy <= y
        model.addConstrs(w[j] <= y[j] for j in range(self.size))
        # lhs
        lhs = [sum(w[j] * self.bi_costs[i][j] for j in range(self.size)) for i in range(self.bi_size)]
        model.addConstrs(lhs[i] >= self.bi_rhss[i] for i in range(self.bi_size))
        model.update()

        # Richard
        model.addConstrs(v[j] * v[j] <= x[j] * y[j] for j in range(self.size))
        # lhs
        lhs_richard = [sum(v[j] * sqrt(self.bi_costs[i][j]) for j in range(self.size)) for i in range(self.bi_size)]
        model.addConstrs(lhs_richard[i] >= sqrt(self.bi_rhss[i]) - EPS_SAFE for i in range(self.bi_size))
        model.update()

        model.optimize()
        tmp_objVal = model.objVal

        model.params.NonConvex = 2
        model.addConstrs(w[j] == x[j] * y[j] for j in range(self.size))
        model.update()

        model.optimize()

        x_sol = [x[j].X for j in range(self.size)]

        y_sol = [y[j].X for j in range(self.size)]

        w_sol = [w[j].X for j in range(self.size)]

        real_lhs = [sum(x_sol[j] * y_sol[j] * self.bi_costs[i][j] for j in range(self.size))
                    for i in range(self.bi_size)]

        gap = [self.bi_rhss[i] - real_lhs[i] for i in range(self.bi_size)]

        if extra_print == 1:
            with np.printoptions(precision=3, suppress=True):
                print("gurobi nonconvex bilinear with Richard:")
                print("obj:")
                print(np.array(model.objVal))
                print("x:")
                print(np.array(x_sol))
                print("y:")
                print(np.array(y_sol))
                print("w:")
                print(np.array(w_sol))
                print("bi_lhs:")
                print(np.array([_.getValue() for _ in lhs]))
                print("gap:")
                print(np.array(gap))

        out_objVal, out_objBound = model.objVal, model.objBound
        del model
        return out_objVal, out_objBound, tmp_objVal

    def print_problem(self):
        with np.printoptions(precision=3, suppress=True):
            print("Problem")
            print("x:")
            print(np.array(self.x_cost))
            print("y:")
            print(np.array(self.y_cost))
            print("b:")
            print(np.array(self.bi_costs))
            print("b-rhs:")
            print(np.array(self.bi_rhss))

    def grb_get_x(self):
        return [self.grb.x[j].X for j in range(self.size)]

    def grb_get_y(self):
        return [self.grb.y[j].X for j in range(self.size)]

    def grb_get_w(self):
        return [self.grb.w[j].X for j in range(self.size)]

    def grb_get_xy(self):
        x_sol = self.grb_get_x()
        y_sol = self.grb_get_y()
        return [x_sol[j] * y_sol[j] for j in range(self.size)]

    def grb_get_bi_gap(self):
        xy_sol = self.grb_get_xy()
        lhs_sol = [sum(xy_sol[j] * self.bi_costs[i][j] for j in range(self.size))
                   for i in range(self.bi_size)]
        gap = [self.bi_rhss[i] - lhs_sol[i] for i in range(self.bi_size)]
        return gap

    def grb_get_val(self):
        return self.grb.m.objVal

    def grb_get_bound(self):
        return self.grb.m.objBound

    def print_grb_sol(self):
        x_sol = [self.grb.x[j].X for j in range(self.size)]

        y_sol = [self.grb.y[j].X for j in range(self.size)]

        real_w = [x_sol[j] * y_sol[j] for j in range(self.size)]

        real_lhs = [sum(real_w[j] * self.bi_costs[i][j] for j in range(self.size))
                    for i in range(self.bi_size)]

        w_sol = [self.grb.w[j].X for j in range(self.size)]

        gap = [self.bi_rhss[i] - real_lhs[i] for i in range(self.bi_size)]

        with np.printoptions(precision=3, suppress=True):
            print("obj:")
            print(np.array(self.grb.m.objVal))
            print("x:")
            print(np.array(x_sol))
            print("y:")
            print(np.array(y_sol))
            print("xy:")
            print(np.array(real_w))
            print("real_lhs:")
            print(np.array(real_lhs))
            print("w:")
            print(np.array(w_sol))
            print("w_lhs:")
            print(np.array([_.getValue() for _ in self.grb.b_lhs]))
            print("gap:")
            print(np.array(gap))

        return self.grb.m.objVal

    # def print_problem_to_file(self, path="./test.npy"):
    #     with open(path, "wb") as file_pointer:
    #         np.save(file_pointer, self.x_cost)
    #         np.save(file_pointer, self.y_cost)
    #         np.save(file_pointer, self.bi_costs)
    #         np.save(file_pointer, self.bi_rhss)
    #
    # def read_problem_from_file(self, path="./test.npy"):
    #     with open(path, "rb") as file_pointer:
    #         x_cost = list(np.load(file_pointer, self.x_cost))
    #         y_cost = list(np.load(file_pointer, self.y_cost))
    #         bi_costs = list(np.load(file_pointer, self.bi_costs))
    #         bi_rhss = list(np.load(file_pointer, self.bi_rhss))
    #     self.add_costs(x_cost, y_cost)
    #     self.add_bilinears(bi_costs, bi_rhss)


class BilinearMosekModel:
    def __init__(self):
        self.m = None
        self.x = None
        self.y = None
        self.w = None  # xy = w
        self.v = None  # xy = v^2
        self.u = None  # min{x,y} = u
        self.b_lhs = None  # bilinear lhs
        self.obj = None

        self.X = None  # [x,y]*[x,y]^T

        self.Z = None  # Extended X with x,y,1
        self.z = None  # 1 within Z

        self._init_mosek_model()

    def _init_mosek_model(self):
        self.m = MF.Model("BPModel")

    def reset_model(self):
        self.__init__()


class BPMosek(BilinearProblem):
    def __init__(self, name=None, OutputFlag=1):
        super().__init__(name, OutputFlag)
        self.msk = None

    def mosek_init(self):
        self.msk = BilinearMosekModel()

    def mosek_reset(self):
        self.msk.reset_model()

    def mosek_init_basic(self):
        self.mosek_init()
        self._mosek_set_var()
        self._mosek_set_var_domain()
        self._mosek_set_obj()
        self._mosek_set_McCormick()
        self._mosek_set_bilinear()

    def _mosek_set_var(self):
        self.msk.m = MF.Model("")
        # self.msk.x = self.msk.m.variable('x', self.size, MF.Domain.inRange(0.0, 1.0))
        # self.msk.y = self.msk.m.variable('y', self.size, MF.Domain.inRange(0.0, 1.0))
        #

        self.msk.Z = self.msk.m.variable('Z', MF.Domain.inPSDCone(2 * self.size + 1))

        self.msk.X = self.msk.Z.slice([0, 0], [2 * self.size, 2 * self.size])
        self.msk.x = self.msk.Z.slice([0, 2 * self.size], [self.size, 2 * self.size + 1])
        self.msk.y = self.msk.Z.slice([self.size, 2 * self.size], [2 * self.size, 2 * self.size + 1])

        self.msk.z = self.msk.Z.slice([2 * self.size, 2 * self.size], [2 * self.size + 1, 2 * self.size + 1])
        self.msk.w = self.msk.Z.pick(range(self.size), range(self.size, 2 * self.size))

    def _mosek_set_var_domain(self):
        self.msk.m.constraint("range", self.msk.Z, MF.Domain.inRange(0.0, 1.0))

    def _mosek_set_obj(self):
        self.msk.obj = MF.Expr.add(MF.Expr.dot(self.x_cost, self.msk.x), MF.Expr.dot(self.y_cost, self.msk.y))
        self.msk.m.objective("obj", MF.ObjectiveSense.Minimize, self.msk.obj)

    def _mosek_set_McCormick(self):
        # tmp_1_vector = [1.0 for _ in range(2 * self.size)]
        tmp_x_y = MF.Expr.vstack(self.msk.x, self.msk.y)
        tmp_v = MF.Expr.repeat(tmp_x_y, 2 * self.size, 1)
        tmp_h = MF.Expr.transpose(tmp_v)
        # X ~ xy <= x
        self.msk.m.constraint("McCormickX", MF.Expr.sub(self.msk.X, tmp_v),
                              MF.Domain.lessThan(0.0))
        # X ~ xy <= y
        self.msk.m.constraint("McCormickY", MF.Expr.sub(self.msk.X, tmp_h),
                              MF.Domain.lessThan(0.0))
        # X ~ xy >= x+y-1
        self.msk.m.constraint("McCormickXY", MF.Expr.add(MF.Expr.sub(self.msk.X, MF.Expr.add(tmp_v, tmp_h)), 1.0),
                              MF.Domain.greaterThan(0.0))
        # X ~ xy >= 0
        # Already set

        self.msk.m.constraint("McCormick1", self.msk.z, MF.Domain.equalsTo(1.0))

    def _mosek_set_bilinear(self):
        for i in range(self.bi_size):
            self.msk.m.constraint("Bilinear"+str(i), MF.Expr.dot(self.bi_costs[i], self.msk.w),
                                  MF.Domain.greaterThan(self.bi_rhss[i]))

    def mosek_opt(self):
        self.msk.m.solve()

    def print_msk_sol(self):
        x_sol = list(self.msk.x.level())
        y_sol = list(self.msk.y.level())
        obj_sol = np.dot(self.x_cost, x_sol) + np.dot(self.y_cost, y_sol)
        real_w = [x_sol[j] * y_sol[j] for j in range(self.size)]
        real_lhs = [sum(real_w[j] * self.bi_costs[i][j] for j in range(self.size))
                    for i in range(self.bi_size)]

        w_sol = list(self.msk.w.level())

        gap = [self.bi_rhss[i] - real_lhs[i] for i in range(self.bi_size)]
        w_lhs = [np.dot(w_sol, self.bi_costs[i]) for i in range(self.bi_size)]

        with np.printoptions(precision=3, suppress=True):
            print("obj:")
            print(np.array(obj_sol))
            print("x:")
            print(np.array(x_sol))
            print("y:")
            print(np.array(y_sol))
            print("xy:")
            print(np.array(real_w))
            print("real_lhs:")
            print(np.array(real_lhs))
            print("w:")
            print(np.array(w_sol))
            print("w_lhs:")
            print(np.array(w_lhs))
            print("gap:")
            print(np.array(gap))

        return obj_sol

if __name__ == "__main__":
    import random
    import types
    import pandas as pd

    n = 10
    itr = 10
    m = 3
    pdList = []
    OutputFlag = 0

    for _ in range(itr):
        x_cost = [random.random() for i in range(n)]
        y_cost = [random.random() for i in range(n)]
        bi_costs = [[random.random() for i in range(n)] for j in range(m)]
        # x_cost = [0.42615123673249633, 0.3346689956407003, 0.05753498671226687, 0.8472632649187982, 0.10667056332148217]
        # y_cost = [0.14636603694915684, 0.3654760714799694, 0.12862825698680647, 0.1542106086938, 0.11667498280195421]
        # bi_cost = [0.23812771801534682, 0.8941770933866809, 0.28623148426147205, 0.3275691346156827, 0.5355427388745678]
        bi_rhss = [sum(bi_costs[j]) / 2 for j in range(m)]

        try:
            BP_test = BPMosek(OutputFlag=OutputFlag)
            BP_test.add_costs(x_cost, y_cost)
            BP_test.add_bilinears(bi_costs, bi_rhss)
            BP_test.print_problem()

            zOpt = BP_test.test_grb_bilinear_opt(OutputFlag=OutputFlag)

            BP_test.grb_init_basic()
            BP_test.grb_opt()
            print("McCormick:")
            zMc = BP_test.print_grb_sol()
            dMc = zMc - zOpt

            BP_test.add_all_minimal_cover_cuts_grb(debug=True)
            BP_test.grb_opt()
            print("Minimal Cover Cut:")
            zMin = BP_test.print_grb_sol()
            dMin = zMin - zOpt

            BP_test.mosek_init_basic()
            BP_test.mosek_opt()
            print("Mosek SDP:")
            zMSK = BP_test.print_msk_sol()
            dMSK = zMSK - zOpt

            if abs(dMc) > 1e-6:
                pdList.append([dMin/dMc, dMSK/dMc])
        except:
            pass

    pdDF = pd.DataFrame(pdList, columns=["RelGap: Minimal Cover", "RelGap: Mosek"])
    pdDF.to_csv("ResultsBilinears.csv", index=False)

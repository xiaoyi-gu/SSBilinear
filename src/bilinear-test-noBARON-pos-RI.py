import bilinear
import warnings
from time import perf_counter
import numpy as np
import random
import pandas as pd
import pathlib
import sys
import pickle


def main():
    target_filename = sys.argv[1]
    target_foldername = str(pathlib.Path(target_filename).parent)
    
    m = int(sys.argv[2])
    n = int(sys.argv[3])
    p = float(sys.argv[4])
    itr = int(sys.argv[5])

    eps_viol = 0.01
    if len(sys.argv) >= 7:
        eps_viol = float(sys.argv[6])

    ind = 0
    if len(sys.argv) >= 8:
        ind = int(sys.argv[7])

    OutputFlag = 0
    eps_heuristic = 0.01
    # eps_zero = 1e-3    # for problem generation
    time_limit = 1800
    # time_limit = 30

    max_search_param = 10
    max_mc_param = max_search_param

    max_mc = int(n * p * max_mc_param)

    if pathlib.Path(target_filename).exists():
        return

    instance_filename = f"instances/{m}-{n}-{p}-{itr}-{eps_viol}-{ind}-pos-RI.pkl"
    fo = None
    fp = None
    x_cost_list, y_cost_list, bi_costs_list, bi_rhss_list = [], [], [], []
    if pathlib.Path(instance_filename).exists():
        fo = open(instance_filename, "rb")
        x_cost_list, y_cost_list, bi_costs_list, bi_rhss_list = pickle.load(fo)

    pdCSV = []
    for _ in range(itr):
        print(f"-- current itr {_} --")
        if fo is None:
            x_cost = [random.random() for i in range(n)]
            y_cost = [random.random() for i in range(n)]
            bi_costs = [[(random.random() if random.random() < p else 0)
                            for i in range(n)] for j in range(m)]
            bi_rhss = [0] * m
            for j in range(m):
                if sum(bi_costs[j]) > 0:
                    bi_rhss[j] = sum(bi_costs[j]) * random.random()
                else:
                    bi_rhss[j] = sum(bi_costs[j]) * (1 + random.random())


            x_cost_list.append(x_cost)
            y_cost_list.append(y_cost)
            bi_costs_list.append(bi_costs)
            bi_rhss_list.append(bi_rhss)
        else:
            x_cost = x_cost_list[_]
            y_cost = y_cost_list[_]
            bi_costs = bi_costs_list[_]
            bi_rhss = bi_rhss_list[_]

        BP_test = bilinear.BPMosek(OutputFlag=OutputFlag)
        BP_test.add_costs(x_cost, y_cost)
        BP_test.add_bilinears(bi_costs, bi_rhss)

        t1 = perf_counter()
        zGRB, zGRB_bound = BP_test.test_grb_bilinear_opt(OutputFlag=OutputFlag, extra_print=1, time_limit=time_limit)
        tGRB = perf_counter() - t1
        print("GRB cpu-time: ", tGRB)

        t1 = perf_counter()
        zRI, zRI_bound, zRI_root = BP_test.test_grb_bilinear_opt_Richard(
            OutputFlag=OutputFlag, extra_print=1, time_limit=time_limit)
        tRI = perf_counter() - t1
        print("RIC cpu-time: ", tRI)

        BP_test.grb_init_basic()
        BP_test.grb_opt()
        print("McCormick:")
        zMc = BP_test.print_grb_sol()

        t1 = perf_counter()
        zPrev = zMc

        n_MC = 0        
        t_grbnlp = 0.0
        t_heuris = 0.0

        for search_itr in range(max_mc):
            xy_sol = BP_test.grb_get_xy()
            bi_gap = BP_test.grb_get_bi_gap()

            t_tmp = perf_counter()
            mc_list = []

            for i in range(BP_test.bi_size):
                if bi_gap[i] > eps_heuristic:
                    label = [None] * n
                    for __ in range(n):
                        if bi_costs[i][__] == 0.0:  # inactive
                            label[__] = -1
                            continue
                        if xy_sol[__] < eps_heuristic:
                            label[__] = 0
                        elif xy_sol[__] > 1 - eps_heuristic:
                            label[__] = 1
                        else:
                            if bi_costs[i][__] > 0:
                                label[__] = 2
                            else:
                                label[__] = int(random.random() < xy_sol[__])    

                    tmp_index_list = [__ for __ in range(n) if label[__] != -1]
                    max_search = max_search_param * len(tmp_index_list)

                    for check_itr in range(max_search):
                        check_mc, tmp_mc = BP_test.check_minimal_cover_generate(label, i)
                        if check_mc == 0:
                            tmp_viol = BP_test._grb_check_valid(tmp_mc)
                            if tmp_viol < -eps_viol:
                                mc_list.append(tmp_mc)
                                n_MC += 1
                            break
                        elif check_mc == -1:
                            tmp_ind_set = [__ for __ in tmp_index_list if label[__] == 1 and bi_costs[i][_] > 0] \
                                        + [__ for __ in tmp_index_list if label[__] == 0 and bi_costs[i][_] < 0]
                            if not tmp_ind_set: break
                            rand_ind = random.choice(tmp_ind_set)
                            if label[rand_ind] == 1:
                                label[rand_ind] = 2
                            else:
                                label[rand_ind] = 1
                        elif check_mc == -2:
                            tmp_ind_set = [__ for __ in tmp_index_list if label[__] == 0 and bi_costs[i][_] > 0] \
                                        + [__ for __ in tmp_index_list if label[__] == 1 and bi_costs[i][_] < 0]
                            if not tmp_ind_set: break
                            rand_ind = random.choice(tmp_ind_set)
                            if label[rand_ind] == 1:
                                label[rand_ind] = 0
                            else:
                                label[rand_ind] = 1
                        elif check_mc == 1:
                            tmp_ind_set = [__ for __ in tmp_index_list if label[__] == 2]
                            if not tmp_ind_set: break
                            tmp_bi_cost = np.array(bi_costs[i])[tmp_ind_set]
                            opt_ind = tmp_ind_set[np.argmin(tmp_bi_cost)]
                            label[opt_ind] = 1
                        else:
                            raise RuntimeError("Unknown return of check_mc from check_minimal_cover()!")

            t_heuris += perf_counter() - t_tmp

            t_tmp = perf_counter()
            for tmp_mc in mc_list:
                BP_test.add_minimal_cover_grb_generated(tmp_mc)
            BP_test.grb_opt()
            zNow = BP_test.print_grb_sol()
            t_grbnlp += perf_counter() - t_tmp

            tMc = perf_counter() - t1
            if ((not mc_list) and zNow - zPrev < 1e-3 * zPrev) or (tMc > time_limit): 
                zPrev = zNow
                break
            zPrev = zNow

        print("Min Cover cpu-time: ", tMc)

        t1 = perf_counter()
        BP_test.test_grb_opt_current(time_limit=time_limit)
        tFin = perf_counter() - t1
        zFin = BP_test.grb_get_val()
        zFin_bound = BP_test.grb_get_bound()

        zOpt = min(zFin, zGRB, zRI)

        pdCSV.append([zMc, zPrev, zGRB, zGRB_bound, zFin, zFin_bound, zOpt,
                      zRI, zRI_bound, zRI_root,
                      n_MC,
                      (zPrev - zMc) / (zOpt - zMc),
                      m, n, p, eps_viol,
                      tGRB, tMc, tFin,
                      tRI,
                      t_grbnlp, t_heuris,
                      ])

    pdDF = pd.DataFrame(pdCSV, columns=["zMc", "zMin-heuristic", "zGRB", "zGRB_bound",
                                        "zFin", "zFin_bound", "zOpt",
                                        "zRI", "zRI_bound", "zRI_root",
                                        "nMC",
                                        "gapClosed-heuristic",
                                        "m", "n", "p", "eps-viol", 
                                        "tGRB", "tMin-heuristic", "tFin",
                                        "tRI",
                                        "tGRBNLP", "tHeuristic",
                                        ])

    pdDF.to_csv(target_filename, index=False)

    if fo is None:        
        fp = open(instance_filename, "wb")
        pickle.dump((x_cost_list, y_cost_list, bi_costs_list, bi_rhss_list), fp)
        fp.close()
    else:
        fo.close()

if __name__ == "__main__":
    main()

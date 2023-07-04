import numpy as np
import pandas as pd

import pathlib
import sys
import os
from itertools import product
import glob

import seaborn as sns
import matplotlib.pyplot as plt


dfList = []
for filename in glob.glob("./results/*-pos-RI.csv"):
    dfList.append(pd.read_csv(filename))

df_orig = pd.concat(dfList)

df_orig['$\\rho_{GRB}$'] = (df_orig['zGRB_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
df_orig['$\\rho_{BC}$'] = (df_orig['zFin_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
df_orig['$\\rho_{MT}$'] = (df_orig['zRI_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100

df_orig['$\\rho_{MTroot}$'] = (df_orig['zRI_root'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100

df_orig['$\\Delta\\rho\\ (\\%)$'] = (df_orig['zFin_bound'] - df_orig['zGRB_bound']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
df_orig['$\\rho_{Heu}\\ (\\%)$'] = df_orig['gapClosed-heuristic'] * 100
# df_orig['$\\rho_{Heu}\\ (\\%)$'] = (df_orig['zMin-heuristic'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100

df_orig.rename(columns={'tMin-heuristic': '$t_{heu}$ (s)', 
                        'tGRB': '$t_{GRB}$',
                        'm': '$m$',
                        'n': '$n$',
                        'p': '$p$',
                        'nMC': '$n_{BC}$',
                    }, inplace=True)
df_orig['$t_{BC}$'] = df_orig['tFin'] + df_orig['$t_{heu}$ (s)']
df_orig['($m$, $p$)'] = tuple(zip(df_orig['$m$'], df_orig['$p$']))

df = df_orig[(df_orig['eps-viol'] == 0.005)]

df['$relgap_{GRB}$'] = (df['zOpt'] - df['zGRB_bound']) / df['zOpt'] * 100
df['$relgap_{BC}$'] = (df['zOpt'] - df['zFin_bound']) / df['zOpt'] * 100



dfsub = df[["($m$, $p$)", "$\\rho_{GRB}$", "$\\rho_{BC}$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('rho-pos-m-p-alter.pdf')


dfsub = df[["($m$, $p$)", "$relgap_{GRB}$", "$relgap_{BC}$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$relgap\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$relgap\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('gap-pos-m-p-alter.pdf')


dfsub = df[["$m$", "$p$", "$\\Delta\\rho\\ (\\%)$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$\\Delta\\rho\\ (\\%)$", x = "$m$")
plt.savefig('rho-pos-p-m.pdf')


dfsub = df[["$m$", "$p$", "$n_{BC}$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$n_{BC}$", x = "$m$")
plt.savefig('n-pos-p-m.pdf')



dfsub = df[["$m$", "$p$", "$\\rho_{Heu}\\ (\\%)$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$\\rho_{Heu}\\ (\\%)$", x = "$m$")
plt.savefig('rhoheu-pos-p-m.pdf')


dfsub = df[["$m$", "$p$", "$t_{heu}$ (s)"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$t_{heu}$ (s)", x = "$m$")
plt.yscale('log')
plt.savefig('t-pos-p-m.pdf')


dfsub = df[["($m$, $p$)", "$\\rho_{MTroot}$", "$\\rho_{Heu}\\ (\\%)$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('rho-pos-m-p-alter-1.pdf')


dfsub = df[["($m$, $p$)", "$\\rho_{GRB}$", "$\\rho_{BC}$", "$\\rho_{MT}$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('rho-pos-m-p-alter-2.pdf')





dfList = []
for filename in glob.glob("./results/*-neg.csv"):
    dfList.append(pd.read_csv(filename))

df_orig = pd.concat(dfList)

df_orig['$\\rho_{GRB}$'] = (df_orig['zGRB_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
df_orig['$\\rho_{BC}$'] = (df_orig['zFin_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
# df_orig['$\\rho_{MT}$'] = (df_orig['zRI_bound'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
# df_orig['$\\rho_{MTroot}$'] = (df_orig['zRI_root'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100

df_orig['$\\Delta\\rho\\ (\\%)$'] = (df_orig['zFin_bound'] - df_orig['zGRB_bound']) / (df_orig['zOpt'] - df_orig['zMc']) * 100
df_orig['$\\rho_{Heu}\\ (\\%)$'] = df_orig['gapClosed-heuristic'] * 100
# df_orig['$\\rho_{Heu}\\ (\\%)$'] = (df_orig['zMin-heuristic'] - df_orig['zMc']) / (df_orig['zOpt'] - df_orig['zMc']) * 100

df_orig.rename(columns={'tMin-heuristic': '$t_{heu}$ (s)', 
                        'tGRB': '$t_{GRB}$',
                        'm': '$m$',
                        'n': '$n$',
                        'p': '$p$',
                        'nMC': '$n_{BC}$',
                    }, inplace=True)
df_orig['$t_{BC}$'] = df_orig['tFin'] + df_orig['$t_{heu}$ (s)']
df_orig['($m$, $p$)'] = tuple(zip(df_orig['$m$'], df_orig['$p$']))

df = df_orig[(df_orig['eps-viol'] == 0.005)]

df['$relgap_{GRB}$'] = (df['zOpt'] - df['zGRB_bound']) / df['zOpt'] * 100
df['$relgap_{BC}$'] = (df['zOpt'] - df['zFin_bound']) / df['zOpt'] * 100



dfsub = df[["($m$, $p$)", "$\\rho_{GRB}$", "$\\rho_{BC}$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('rho-neg-m-p-alter.pdf')


dfsub = df[["($m$, $p$)", "$relgap_{GRB}$", "$relgap_{BC}$"]]
dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$relgap\\ (\\%)$", var_name = " ")
plt.figure(figsize=(10,5))
g = sns.barplot(data=dfsub, hue = " ", y = "$relgap\\ (\\%)$", x = "($m$, $p$)")
g.get_legend().set_title(None)
plt.savefig('gap-neg-m-p-alter.pdf')


dfsub = df[["$m$", "$p$", "$\\Delta\\rho\\ (\\%)$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$\\Delta\\rho\\ (\\%)$", x = "$m$")
plt.savefig('rho-neg-p-m.pdf')


dfsub = df[["$m$", "$p$", "$n_{BC}$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$n_{BC}$", x = "$m$")
plt.savefig('n-neg-p-m.pdf')



dfsub = df[["$m$", "$p$", "$\\rho_{Heu}\\ (\\%)$"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$\\rho_{Heu}\\ (\\%)$", x = "$m$")
plt.savefig('rhoheu-neg-p-m.pdf')


dfsub = df[["$m$", "$p$", "$t_{heu}$ (s)"]]
plt.figure(figsize=(8,5))
sns.barplot(data=dfsub, hue = "$p$", y = "$t_{heu}$ (s)", x = "$m$")
plt.yscale('log')
plt.savefig('t-neg-p-m.pdf')


# dfsub = df[["($m$, $p$)", "$\\rho_{MTroot}$", "$\\rho_{Heu}\\ (\\%)$"]]
# dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
# plt.figure(figsize=(10,5))
# g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
# g.get_legend().set_title(None)
# plt.savefig('rho-neg-m-p-alter-1.pdf')


# dfsub = df[["($m$, $p$)", "$\\rho_{GRB}$", "$\\rho_{BC}$", "$\\rho_{MT}$"]]
# dfsub = dfsub.melt(id_vars = "($m$, $p$)", value_name = "$\\rho\\ (\\%)$", var_name = " ")
# plt.figure(figsize=(10,5))
# g = sns.barplot(data=dfsub, hue = " ", y = "$\\rho\\ (\\%)$", x = "($m$, $p$)")
# g.get_legend().set_title(None)
# plt.savefig('rho-neg-m-p-alter-2.pdf')
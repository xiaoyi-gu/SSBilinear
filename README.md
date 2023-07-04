# SSBilinear - Solving *S*parse *S*eparable *Bilinear* Programs Using Lifted Bilinear Cover Inequalities
Computational codes for Gu X, Dey SS, Richard JPP (2023) "Solving sparse separable bilinear programs using lifted bilinear cover inequalities".

## Prerequisites (*Settings in Paper*)
* Python 3 (*Python 3.9.7*)
* Gurobi (*Gurobi 9.5.0*)

## Latest Release and Results 
[Link](https://github.com/xiaoyi-gu/SSBilinear/releases/latest)

## Installation
For Local Execution:
* Download the codes and instances, and extract into a proper folder.
Warning: local and unparallelized execution will take an extended period of time.

For HTCondor Execution:
* Download the codes and instances, and extract into a proper folder.
* (Optional: setup virtual environment if needed): edit `mosek-env-setup.sh` to fit your condor settings and run `bash mosek-env-setup.sh`.
* Modify several `make_*.sh` to fit your condor settings and/or gurobi settings.
* Modify several `htcondor_*.cmd` to fit your condor settings and/or python settings.

## Execution
For Local Execution:
* Execute
```bash
python python_run.py
```

For HTCondor Execution:
* Execute
```bash
bash condor_run.sh
```
* Alternatively, execute `condor_submit htcondor_*.cmd` for the proper `htcondor_*.cmd`.

## Figures
For Figure Generation:
* Collect results from execution. Alternatively, use the results we generated
* Execute `python generate_figures.py`. This can be doen within a separate python environment, as long as proper packages like `seaborn` are installed.

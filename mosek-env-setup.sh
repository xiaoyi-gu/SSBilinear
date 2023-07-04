#!/bin/bash
# /home/xgu74/.julia/conda/3/bin/python3 -m venv mosek-env
python3 -m venv mosek-env
source ./mosek-env/bin/activate
pip3 install --upgrade pip
pip3 install pandas
pip3 install numpy
pip3 install xlsxwriter
pip3 install xlrd
pip3 install mosek
pip3 install pyomo
python -m pip install gurobipy==9.5.0
deactivate
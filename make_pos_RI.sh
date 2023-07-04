#!/bin/bash
source mosek-env/bin/activate
GUROBI_HOME="/opt/gurobi/gurobi950/" GRB_LICENSE_FILE="/opt/gurobi/gurobi950/linux64/gurobi-client.lic" python src/bilinear-test-noBARON-pos-RI.py $1 $2 $3 $4 $5 $6 $7 > $1.log 2>&1
deactivate
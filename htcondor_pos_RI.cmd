universe = vanilla
getenv = true
should_transfer_files = IF_NEEDED
executable = make_pos_RI.sh
arguments = results/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos-RI.csv $(M) $(N) $(P) $(I) $(VIOL) $(IND)
Log = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos-RI.csv.log
output = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos-RI.csv.out
error = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos-RI.csv.error
notification = error
notification = complete
request_memory = 16GB
request_cpus = 1
queue M,N,P,I,VIOL,IND from qList.txt
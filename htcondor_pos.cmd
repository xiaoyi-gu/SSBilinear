universe = vanilla
getenv = true
should_transfer_files = IF_NEEDED
executable = make_pos.sh
arguments = results/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos.csv $(M) $(N) $(P) $(I) $(VIOL) $(IND)
Log = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos.csv.log
output = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos.csv.out
error = ./HTCondorLog/$(M)-$(N)-$(P)-$(I)-$(VIOL)-$(IND)-pos.csv.error
notification = error
notification = complete
request_memory = 16GB
request_cpus = 1
queue M,N,P,I,VIOL,IND from qList.txt
from subprocess import *

# Run first run with no filters
proc1 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -sp -pq 1 -mcr 0 -op Bueno/no_filter', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc1.communicate()
#Get plots for run with no filters
proc11 = Popen('python withinCancerEvaluation.py -cf Bueno/BUENO_config_deep.json -op Bueno/no_filter', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc11.communicate()

# Run first run with only permuted q-value
proc2 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -pq 0.1 -mcr 0 -op Bueno/pq', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc2.communicate()
#Get plots for run with only permuted q-value
proc21 = Popen('python withinCancerEvaluation.py -cf Bueno/BUENO_config_deep.json -op Bueno/pq', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc11.communicate()

# Run first run with only shallow coincidence filter
proc3 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -pq 1 -mcr 0.001 -op Bueno/mcr', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc3.communicate()
#Get plots for run with only shallow coincidence filter
proc31 = Popen('python withinCancerEvaluation.py -cf Bueno/BUENO_config_deep.json -op Bueno/mcr', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc11.communicate()

# Run first run with both filters
proc4 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -sp -pq 0.1 -mcr 0.001 -op Bueno/pq_mcr', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc4.communicate()
#Get plots for run with both filters
proc41 = Popen('python withinCancerEvaluation.py -cf Bueno/BUENO_config_deep.json -op Bueno/pq_mcr', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc11.communicate()

# Get Plots From Across all Runs
proc5 = Popen('python crossCancerEvaluation.py -cf Bueno/BUENO_config_deep.json -s Bueno/Summary_Paths.xlsx', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
out = proc11.communicate()



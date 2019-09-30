# For running the OncoMerge program generally
from subprocess import *
import pandas as pd
import os
import csv
from multiprocessing import Pool, cpu_count
import oncoMerge_utils as omu
import time
import json
from pathlib import Path


### Enter Parameters Here

# Unique Name for the current run. Reccomend name/abbreviation of tumor type. Must be passed as a list of strings. Even when there is only one element.
Run_Names = sorted(['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC'])
# Minimun Mutation Frequency (recommend 0.05):
mutation_frequency =  0.05

# p-value to use for Permutation Analysis (reccommend 0.1):
#perm_pv = 0.1
# Vary permpv and coinc
'''mc_pv_combos = [(0.001,0.1),(0,0.1),(0.001,1),(0,1)]
Run_Names1 = []
for mc_pv_combo in mc_pv_combos:
    mc, perm_pv = mc_pv_combo
    #choose where configs go
    for Run_Name in Run_Names:
        Run_Name1 = Run_Name + '_minCoinc_'+str(mc)+'_permPv_'+str(perm_pv)
        Run_Names1.append(Run_Name1)


#Build Configs
#Creates path strings for config to fill params later.
for Run_Name1 in Run_Names1:
    Run_Name, str1, mc, str2, perm_pv = Run_Name1.split('_')
'''
config_paths = {}
for Run_Name in Run_Names:
    config_paths[Run_Name] = 'configs/final/'+Run_Name+'_config.JSON'

min_loci_genes = 2
mc = 0.001
perm_pv = 0.1
for Run_Name in Run_Names:
    Run_Name1 = Run_Name
    mc = float(mc)
    perm_pv = float(perm_pv)
    params = {}
    if Run_Name == 'SKCM':
        p = 'TM'
    else:
        p = 'TP'
    params['conversion_file_path'] = 'OncoMerge_input_g2e_converter.csv'
    params['gistic_path'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/'
    #params['Dels'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt'
#Amplifications File:
    #params['Amps'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt'
# Gene Data File:
    #params['gene_data_path'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt'
# Title of Column containing entrez ids in the gistic Gene Data file (all_data_by_genes.txt) should match title of entrez id column in symbol_conversion_path:
    params['label_name'] = 'Locus ID'
#all_Thresholded File:
    #params['thresh'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt'
# Somatic Mutations File:
    params['som_mut_path'] = 'mutations_mc3/'+Run_Name+'_somMutMC3.csv'
# MutSig2CV File:
    params['mutsig2CV_path'] = 'sig2cv/'+Run_Name+'_sig2cv.csv'
# Where to store output:
    params['output_path'] = 'output/final/oncoMerged_'+Run_Name1
    params['log'] = 'output/final/log/' +Run_Name1+ '.log'
    params['CNALoci'] = params['output_path']+ '/oncoMerge_CNA_loci_.csv'
    params['oncoMerged'] = params['output_path']+ '/oncoMerge_mergedMuts.csv'
    params['LofActSig'] = params['output_path']+ '/oncoMerge_ActLofPermPV.csv'
    params['config_path'] = config_paths[Run_Name1]
    params['min_mut_freq'] = mutation_frequency
    params['perm_pv'] = perm_pv
    params['min_loci_genes'] = min_loci_genes
    params['min_coinc'] = mc
    if not os.path.exists(params['config_path']):
        Path(params['config_path']).touch()
    for i in range(sum('/' == i for i in params['output_path'])):
        if not os.path.exists('/'.join(params['output_path'].split('/')[:i+2])):
            os.mkdir('/'.join(params['output_path'].split('/')[:i+2]))                   #
    with open(params['config_path'], 'w') as cf:
        json.dump(params, cf)

# Run OncoMerge
def RDProc(a):
    start = time.time()
    with open(a[-1],'r') as cf:
        params = json.loads(cf.read())
    print(' '.join(a))
    logFile = open(params['log'],'w')
    RDPproc = Popen(' '.join(a), stdout=logFile, stderr=STDOUT, shell = True)
    output = RDPproc.communicate()
    end = time.time()
    print(end-start)



def main(Run_Names1):
    runMe = []
    for Run_Name1 in Run_Names1:
        with open(config_paths[Run_Name1],'r') as cf:
            params = json.loads(cf.read())
        if not os.path.exists(params['oncoMerged']):
           #runMe.append(['python oncoMerge_9_13.py', ' -cf', config_paths[Run_Name1]])
           runMe.append(['python oncoMerge.py', ' -cf', config_paths[Run_Name1]])
    cpus = 2
    # cpus = cpu_count()
    print('There are %d CPUs available.' % cpus)
    print(runMe)
    pool = Pool(processes=cpus)
    #return runMe
    pool.map(RDProc, runMe)
    pool.close()
    pool.join()
    return runMe
#for a in runMe:
#    RDProc(a)

if __name__ == "__main__":
        #print('Creating conversion File')
        #logFile = open('output/final/log/conversionConstrution.log','w')
        #CC = Popen('python make_AmpDel_Genelist.py', stdout=logFile, stderr=STDOUT, shell = True)
        #CC.wait()
        #print('Done with conversion File!')
        runMe = main(Run_Names)
        #
        print('Starting Tests')
        for Run_Name in Run_Names:
            logFile = open('output/final/log/Test_'+Run_Name+'.log','w')
            oMtest =  Popen('python test_oncoMerge.py -p1 '+'mcpv_output/oncoMerged_'+Run_Name+'_minCoinc_0.001_permPv_0.1/oncoMerge_mergedMuts.csv'+' -p2 '+'output/final/oncoMerged_'+Run_Name+'/oncoMerge_mergedMuts.csv'+' -rn '+ Run_Name, stdout=logFile, stderr=STDOUT, shell = True)
            oMtest.wait()
        print('Done with Filter Summary Table')
        #
        print('Starting Summary Tables')
        logFile = open('output/final/log/BuildSummaryTable.log','w')
        Assess_Value_Added = Popen('python Build_OM_Summary_final.py -o output/final/', stdout=logFile, stderr=STDOUT, shell = True)
        Assess_Value_Added.wait()
        print('Done with om Summary!')
        '''print('Starting Filter Summary Table')
        logFile = open('output/final/log/oM_Summary_Table.log','w')
        oMsum = Popen('python Build_oncoMerge_SummaryTable_8_15.py', stdout=logFile, stderr=STDOUT, shell = True)
        oMsum.wait()
        print('Done with Filter Summary Table')
        '''
        print('Starting Correlation Analysis')
        logFile = open('output/final/log/Corellation_analysis.log','w')
        CorAn = Popen('python output/Correlation_Analysis.py -o output/final/', stdout=logFile, stderr=STDOUT, shell = True)
        CorAn.wait()
        print('Done with Correlation Analysis')
        '''
        print('Starting True Neg Calc')
        logFile = open('output/final/log/TrueNegs.log','w')
        TrueNegs = Popen('python Calc_True_Neg.py', stdout=logFile, stderr=STDOUT, shell = True)
        TrueNegs.wait()
        print('Done with True Negs')
        '''

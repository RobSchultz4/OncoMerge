# Compare and Contrast Mutational Landscapes of different oncoMerge Runs

# Preamble
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from matplotlib import pyplot as plt
import argparse
import numpy as np
import json


# Read in Command Line Options   # need to choose how to lead in the necessary data
parser = argparse.ArgumentParser(description='Used to choose a cancer type')
parser.add_argument('-s', '--sum_path', help='Path to a csv file containing 2 columns: the first is the run name and the second is the path to the oncoMerge Summary Matrix for the corresponcing run name.', type = str)
parser.add_argument('-cf', '--config_path', help='Path to config file used for oncoMerge', type = str)
args = parser.parse_args()

params = args.__dict__
if args.config_path:
    with open(args.config_path, "r") as cfg:
        tmp = json.loads(cfg.read())
        for i in tmp:
            params[i] = tmp[i]

# Aggregate oncoMerge Summary Matrices
print('Aggregate Summary Tables')
sumTablePaths = {}
with open(args.sum_path,'r') as inFile:
    lines = [line for line in inFile.read().split('\n') if not line =='']
    for line in lines:
        splitUp = line.split(',')
        sumTablePaths[splitUp[0]] = splitUp[1]

sumTables = {}
for runName in sumTablePaths.keys():
    table1 = pd.read_csv(sumTablePaths[runName], index_col = 0)
    sumTables[runName]  = table1.loc[~pd.isnull(table1['Final_mutation_type'])]

summaryMatrix = pd.concat(sumTables)
tmpcol = summaryMatrix.reset_index().ix[:,0].to_frame()
tmpcol.columns = ['Tumor_type']
tmpcol.index = summaryMatrix.index
summaryMatrix = pd.concat([summaryMatrix,tmpcol],axis = 1)
summaryMatrix.index.levels[0].name = 'Run_name'
summaryMatrix.index.levels[1].name = 'Entrez'
summaryMatrix['Final_freq'] = summaryMatrix['Final_freq'].astype(float)
summaryMatrix['Delta_over_PAM'] = summaryMatrix['Delta_over_PAM'].replace(np.nan, 0)
summaryMatrix['CNA_freq'] = summaryMatrix['CNA_freq'].replace(np.nan, 0)
summaryMatrix['PAM_freq'] = summaryMatrix['PAM_freq'].replace(np.nan, 0)
summaryMatrix['PAM_only'] = summaryMatrix['Final_freq'] - summaryMatrix['CNA_freq']
summaryMatrix['Coincidence_rate'] = summaryMatrix['PAM_freq'] + summaryMatrix['CNA_freq'] - summaryMatrix['Final_freq']
summaryMatrix['CNA_only'] = summaryMatrix['Final_freq'] - summaryMatrix['PAM_freq']

# Cross Cancer Evaluations

# Set Symbols as index
summaryMatrix21 = summaryMatrix.reset_index().set_index('Symbol')

# What genes are most common in cancers?
# Get the number of cancers each gene is important in

#crossCancer = pd.DataFrame(summaryMatrix['Symbol'].value_counts())
crossCancer = pd.DataFrame(summaryMatrix['Symbol'].value_counts())
crossCancer.columns = ['Num_tumor_types']
# Add list of cancers
for mut in list(set(summaryMatrix['Symbol'])):
    if crossCancer.loc[mut,'Num_tumor_types'] > 1:
        crossCancer.at[mut,'Tumor_types'] = ', '.join(set(summaryMatrix21.loc[mut,'Tumor_type']))
    else:
        crossCancer.at[mut,'Tumor_types'] = summaryMatrix21.loc[mut,'Tumor_type']

# Add Mutypes
crossCancer['Final_mutation_types'] = ''
for gene in crossCancer.index:
    if crossCancer.loc[gene,'Num_tumor_types' ] > 1:
        mTypes1 = ', '.join(list(summaryMatrix21.loc[gene,'Final_mutation_type']))
    else:
        mTypes1 = summaryMatrix21.loc[gene,'Final_mutation_type']
    crossCancer.at[gene, 'Final_mutation_types'] = mTypes1


crossCancer.to_csv(params['output_path']+'/Muts_Across_TumorTypes.csv')

# How many genes do each cancer type share with eachother?
TTtable = pd.DataFrame(0,index = sumTablePaths.keys(), columns = sumTablePaths.keys())

with PdfPages('TumorTypeClusters_sharedGenes.pdf') as pdf
    ax = sns.clustermap(TTtable)
    plt.title('Tumor Type Clusters')

for runName1 in sumTablePaths.keys():
    for runName2 in sumTablePaths.keys():
        for l1 in crossCancer.loc[crossCancer['Num_tumor_types']>1]['Tumor_types']:
            if runName1 in l1:
                if runName2 in l1:
                    TTtable.loc[runName1,runName2] +=1

# What Cancers are paired together the most often?
ttsims = pd.DataFrame(index = sumTablePaths.keys(), columns = ['Most_similar_tumor_type', 'Genes'])

TTtable.to_csv(params['output_path']+'/TumorType_Similarity.csv')
for runName in sumTablePaths.keys():
    paired = list(TTtable.loc[runName][TTtable.loc[runName] == TTtable.loc[runName,[i for i in sumTablePaths.keys() if not i==runName]].max()].index)
    ttsims.at[runName,'Most_similar_tumor_type'] = paired
    ttsims.at[runName,'Genes'] = []
    for pair in paired:
        pairdf1 = summaryMatrix.loc[[runName,pair]]
        ttsims.at[runName,'Genes'].append(list(set([i for i in pairdf1['Symbol'] if pairdf1['Symbol'].value_counts().loc[i] ==2 ])))


ttsims.to_csv(params['output_path']+'/TTsim_summary.csv')

crossCancer.loc[[i for i in crossCancer.index if any(mt in crossCancer.loc[i,'Final_mutation_types'] for mt in ['Act','CNAamp'])  and any(mt in crossCancer.loc[i,'Final_mutation_types'] for mt in ['LoF', 'CNAdel'])]].to_csv(params['output_path']+'/Context_dependent_mutations.csv')



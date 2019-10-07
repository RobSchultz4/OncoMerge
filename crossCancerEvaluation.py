# Generate Plots comparing oncoMerge output across cancers

# Preamble
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from matplotlib import pyplot as plt
#from subprocess import *
#import os
#from scipy.stats import hypergeom as hyg
#import numpy as np
#from scipy.stats import pearsonr
#from itertools import combinations
#from math import log10
import argparse
#import json
#import dill
#from copy import deepcopy
#import pickle as pkl


# Read in Command Line Options   # need to choose how to lead in the necessary data
parser = argparse.ArgumentParser(description='Used to choose a cancer type')
parser.add_argument('-s', '--sum_path', help='Path to a csv file containing 2 columns: the first is the run name and the second is the path to the oncoMerge Summary Matrix for the corresponcing run name.', type = str)
parser.add_argument('-cf', '--config_path', help='Path to config file used for oncoMerge', type = str)
args = parser.parse_args()

# Define Functions


# Define Parameters
MutationTypes = ['PAM','CNAamp','Act','CNAdel','LoF']
mf = 0.05

# Aggregate oncoMerge Summary Matrices
with open(args.sum_path,'r') as inFile:
    sumTablePaths = {i.split(',')[0]:i.split(',')[1] for i in inFile.read().split('\n')}

for runName in sumTablePaths.keys():
    table1 = pd.read_csv(sumTablePaths[runName])
    sumTables[runName]  = table1.loc[~pd.isnull(table1['Final_mutation_type'])]

summaryMatrix = pd.concat(sumTables)
tmpcol = summaryMatrix.reset_index().ix[:,0].to_frame()
tmpcol.columns = ['Tumor_type']
tmpcol.index = summaryMatrix.index
summaryMatrix = pd.concat([summaryMatrix,tmpcol],axis = 1)
summaryMatrix['Final_freq'] = summaryMatrix['Final_freq'].astype(float)
summaryMatrix['Delta_over_PAM'] = summaryMatrix['Delta_over_PAM'].replace(np.nan, 0)
summaryMatrix['CNA_freq'] = summaryMatrix['CNA_freq'].replace(np.nan, 0)
summaryMatrix['PAM_freq'] = summaryMatrix['PAM_freq'].replace(np.nan, 0)
summaryMatrix['PAM_only'] = summaryMatrix['Final_freq'] - summaryMatrix['CNA_freq']
summaryMatrix['Coincidence_rate'] = summaryMatrix['PAM_freq'] + summaryMatrix['CNA_freq'] - summaryMatrix['Final_freq']
summaryMatrix['CNA_only'] = summaryMatrix['Final_freq'] - summaryMatrix['PAM_freq']


#################
## Begin Plots ##
#################

### Make Bar plot and raw csv of number of each CNAtype in OncoMerged file
mutTypes = {}
#maxmuts = 400
ind = range(len(MutationTypes))
with PdfPages(args.output_path+ '/All_MutTypes_BarPlot.pdf') as pdf:
    for runName in sumTables.keys():                                                                #
        mutTypes[runName] = sumTables[runName].groupby('Final_mutation_type').count()['Final_freq']
        for mt in MutationTypes:
            if not mt in mutTypes[runName].index:
                mutTypes[runName].loc[mt] = 0
        mutTypes[runName] = mutTypes[runName].loc[['PAM','Act','CNAamp','LoF','CNAdel']]
        plt.bar(ind,mutTypes[runName], color = ['r','b','g','c','m'])
        plt.xticks(ind, mutTypes[runName].index)
        plt.title(runName+', Total Genes: '+str(len(summaryMatrix.loc[runName])))
        pdf.savefig()
        plt.close()


mutTypeDf = pd.concat(mutTypes).unstack()
mutTypeDf.to_csv(args.output_path+ '/All_MutTypes_BarPlot_raw.csv')

# Frequency Violin
# Plot Violin plots of the frequency of mutation types for each cancer

with PdfPages(args.output_path+ '/All_FreqViolins_byTumorType_point.pdf') as pdf:
    for mt in MutationTypes:
        ax = sns.violinplot(x = 'Tumor_type', y = 'Final_freq', scale = 'count',inner = 'stick', order = list(mutTypeDf.sort_values('PAM').index) , data = summaryMatrix.loc[summaryMatrix['Final_mutation_type']==mt,['Final_freq','Tumor_type','Final_mutation_type']],axis = 1)
        #ax = sns.swarmplot( x = 'Cancers', y = 'Mutation Frequency', color = 'k', size = 3, order = list(dv.sort_values('PAM').index) , data=smsallthings[sms[i]])#,axis = 1)
        plt.title(mt)
        plt.ylim(float(mf), 1)
        plt.xticks(rotation=90)
        pdf.savefig()
        plt.close()

#Plots Violins for each Cancer type showing distribution of mutation frequency for different Mutation Types.
with PdfPages(args.output_path+ '/All_FreqViolins_byMutType_swarm.pdf') as pdf:
    for tt in summaryMatrix.index.levels[0]:
        data = summaryMatrix.loc[summaryMatrix['Tumor_type'] == tt, ['Final_freq','Final_mutation_type']]
        ax = sns.violinplot(x = 'Final_mutation_type', y = 'Final_freq', scale = 'count', inner = None, order = ['PAM','Act','CNAamp','LoF','CNAdel'], data = data)
        ax = sns.swarmplot( x = 'Final_mutation_type', y = 'Final_freq', color = 'k',     size = 2,     order = ['PAM','Act','CNAamp','LoF','CNAdel'], data = data) #, ax = ax)
        plt.ylim(float(mf), 1)
        plt.title('Frequency of Mutation in ' + tt)
        pdf.savefig()
        plt.close()

# Make Freq Extension Violin plots for each cancer

freqEx_order = pd.Series(index = summaryMatrix.index.levels[0] )
for tt in summaryMatrix.index.levels[0]:
    freqEx_order.loc[tt] = summaryMatrix.loc[summaryMatrix['Tumor_type']==tt]['Delta_over_PAM'].mean()

fe_or = ['All'] + list(freqEx_order.sort_values(ascending = False).index)
#fe_or = list(freqEx_order.sort_values(ascending = False).index)
s1 = summaryMatrix.reset_index()
s1.ix[:,0] = 'All'
s1 = s1.set_index(['level_0','level_1'])
data1 = pd.concat([summaryMatrix,s1])
with PdfPages(args.output_path+ '/All_FreqExtens_Violins_tt.pdf') as pdf:
    ax = sns.violinplot(x = 'Tumor_type', y = 'Delta_over_PAM', data = data1, cut = 0, inner = None, scale = 'count', order = fe_or)
    plt.ylim(0, 1)
    #plt.title('oncoMerge Value Added')
    plt.xticks(rotation=90)
    pdf.savefig()
    plt.close()

recoveredMuts_df = summaryMatrix.loc[(summaryMatrix['PAM_freq'] < 0.05) & (summaryMatrix['CNA_freq'] < 0.05)]
recovered_MutTypes_Raw_df = pd.DataFrame(columns = ['Act','LoF'])
for tt in summaryMatrix.index.levels[0]:
    for sm in ['Act', 'LoF']:
        if any(tt in i for i in recoveredMuts_df.index):
            recovered_MutTypes_Raw_df.loc[tt,sm] = sum(recoveredMuts_df.loc[tt]['Mutation Type'] == sm)
        else:
            recovered_MutTypes_Raw_df.loc[tt,sm] = 0

recovered_MutTypes_Raw_df.to_csv(args.output_path + '/All_recovered_MutTypes_Raw.csv')

# Create Stacked bar plot of the CNAs, PAMs, and coincidence for all tts
with PdfPages(args.output_path+'/All_StackedBars.pdf') as pdf:
    for tt in summaryMatrix.index.levels[0]:
        df1 = summaryMatrix.loc[tt]
        df1 = df1.loc[list(df1['Final_mutation_type'].sort_values().index)]
        N = len(df1)
        ind = np.arange(N)
        p1 = plt.bar(ind,df1['PAM_only'], color = 'y')
        p2 = plt.bar(ind,df1['Coincidence_rate'], bottom = df1['PAM_only'] , color = 'g')
        p3 = plt.bar(ind,df1['CNA_only'], bottom =df1['Coincidence_rate']+ df1['PAM_only'], color = 'b')
        p4 = plt.plot([-1,len(ind)], [0.05,0.05], 'r', label='Mutation Frequency Threshold', linewidth = 1/len(ind))
        plt.ylabel('Frequency')
        plt.title('oncoMerged Mutation Frequencies in '+tt)
        plt.xticks(ind,tuple(df1['Symbol']), rotation = 90)
        plt.legend(tuple([p1[0],p2[0],p3[0],p4[0]]), tuple(['PAM only', 'Both','CNA only', 'Threshold']))
        plt.tight_layout()
        pdf.savefig()
        plt.close()












# Cross Cancer Evaluations

cancers = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

# Currently assumes summaryMatrix is loaded up
summaryMatrix21 = summaryMatrix.reset_index().set_index('Symbol')

# What genes are most common in cancers?
# Get the number of cancers each gene is important in
crossCancer = pd.DataFrame(summaryMatrix['Symbol'].value_counts())
crossCancer.columns = ['Num Tumor Types']
# Add list of cancers
for mut in list(set(summaryMatrix['Symbol'])):
    if crossCancer.loc[mut,'Num Tumor Types'] > 1:
        crossCancer.at[mut,'Cancers'] = ', '.join(set(summaryMatrix21.loc[mut,'Tumor Type']))
    else:
        crossCancer.at[mut,'Cancers'] = summaryMatrix21.loc[mut,'Tumor Type']

# Add Mutypes
crossCancer['Mutation Types'] = ''
for gene in crossCancer.index:
    if crossCancer.loc[gene,'Num Tumor Types' ] > 1:
        mTypes1 = ', '.join(list(summaryMatrix21.loc[gene,'Mutation Type']))
    else:
        mTypes1 = summaryMatrix21.loc[gene,'Mutation Type']
    crossCancer.at[gene, 'Mutation Types'] = mTypes1


crossCancer.to_csv('output/final_10_1/analysis/Muts_Across_TumorTypes.csv')

# How many genes do each cancer type share with eachother?
TTtable = pd.DataFrame(0,index = cancers, columns = cancers)

for cancer1 in cancers:
    for cancer2 in cancers:
        for l1 in crossCancer.loc[crossCancer['Num Tumor Types']>1]['Cancers']:
            if cancer1 in l1:
                if cancer2 in l1:
                    TTtable.loc[cancer1,cancer2] +=1

# What Cancers are paired together the most often?
ttsims = pd.DataFrame(index = cancers, columns = ['Most Similar Tumor Type', 'Genes'])

TTtable.to_csv('output/final_10_1/analysis/TumorType_Similarity.csv')
for cancer in cancers:
    paired = list(TTtable.loc[cancer][TTtable.loc[cancer] == TTtable.loc[cancer,[i for i in cancers if not i==cancer]].max()].index)
    ttsims.at[cancer,'Most Similar Tumor Type'] = paired
    ttsims.at[cancer,'Genes'] = []
    for pair in paired:
        pairdf1 = summaryMatrix.loc[[cancer,pair]]
        ttsims.at[cancer,'Genes'].append(list(set([i for i in pairdf1['Symbol'] if pairdf1['Symbol'].value_counts().loc[i] ==2 ])))


ttsims.to_csv('output/final_10_1/analysis/TTsim_summary.csv')

crossCancer.loc[[i for i in crossCancer.index if any(mt in crossCancer.loc[i,'Mutation Types'] for mt in ['Act','CNAamp'])  and any(mt in crossCancer.loc[i,'Mutation Types'] for mt in ['LoF', 'CNAdel'])]].to_csv('output/final_10_1/analysis/Context_dependent_mutations.csv')



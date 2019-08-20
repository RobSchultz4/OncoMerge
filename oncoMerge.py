##########################################################
## OncoMerge:  oncoMerge.py                             ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: https://github.com/plaisier-lab/OncoMerge   ##
## @Author:  Chris Plaisier, Rob Schultz                ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

#####################
## Import packages ##
#####################
import json
import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import itertools
import sys
import os

##################################
## Read in command line options ##
##################################
parser = argparse.ArgumentParser(description='Used to choose a cancer type')
parser.add_argument('-cf', '--config_path', help='Path to JSON encoded configuration file, overrides command line parameters', type = str)
parser.add_argument('-gp', '--gistic_path', help='Path to GISTIC output folder', type = str)
parser.add_argument('-dp', '--del_path', help='Path to GISTIC deletion file (default = del_genes.conf_99.txt)', default = 'del_genes.conf_99.txt', type = str)
parser.add_argument('-ap', '--amp_path', help='Path to the GISTIC amplification file (default = amp_genes.conf_99.txt)', default = 'amp_genes.conf_99.txt', type = str)
parser.add_argument('-gdp', '--gene_data_path', help='Path to the GISTIC gene data file (default = all_data_by_genes.txt)', default = 'all_data_by_genes.txt', type = str)
parser.add_argument('-cfp', '--conversion_file_path', help='Supply hand annotated file to convert gene symbols to Entrez IDs (default does not import hand annotated conversion file and instead uses conversion embedded in GISTIC output files).', type = str)
parser.add_argument('-ln', '--label_name', help='Label for Entrez ID column in GISTIC gene data file (default = \'Gene ID\')', type = str, default='Gene ID')
parser.add_argument('-tp', '--thresh_path', help='Path to the GISTIC all_thresholded file (all_thresholded.by_genes.txt)', default = 'all_thresholded.by_genes.txt', type = str)
parser.add_argument('-smp', '--som_mut_path', help='Path to the somatic mutations file (CSV matrix where columns are patients and genes are rows) [0 = not mutated, and 1 = mutated]', type = str)
parser.add_argument('-mscv', '--mutsig2_cv_path', help='Path to a MutSig2CV output file', type = str)
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default='.')
parser.add_argument('-mmf', '--min_mut_freq', help='Minimum frequency of mutation (range = 0-1; default = 0.05)', type = float, default='0.05')
parser.add_argument('-pp', '--perm_pv', help='Permuted p-value FDR BH corrected cutoff (default = 0.1)', type = float, default='0.1')
args = parser.parse_args()

#######################
## Define parameters ##
#######################
params = args.__dict__
if args.config_path:
    with open(args.config_path, "r") as cfg:
        tmp = json.loads(cfg.read())
        for i in tmp:
            params[i] = tmp[i]
if (not params['gistic_path']) or (not params['som_mut_path']) or (not params['mutsig2CV_path']):
    parser.print_help()
    sys.exit(1)

# Set minimum coincidence rate
params['min_coinc'] = 0.001

##################
## Load up data ##
##################
# Create conversion series for gene symbol to Entrez ID
if not params['conversion_file_path']:
    n1 = pd.read_csv(params['gistic_path']+'/'+params['gene_data_path'],index_col=0,sep='\t',usecols=[0,1])
    n1.index = [i.split('|')[0] for i in n1.index]
    n1 = n1.drop_duplicates()
    n1 = n1[params['label_name']]
else:
    n1 = pd.read_csv(params['conversion_file_path'],index_col=0)[label_name].apply(int)

# load up significantly mutated genes
mutSig2CV = pd.read_csv(params['mutsig2CV_path'],index_col=1)
mutSig2CV = mutSig2CV.loc[mutSig2CV.index.map(lambda x: x in n1.index)]
mutSig2CV.index = mutSig2CV.index.map(lambda x: n1.loc[x])
mutSig2CV = mutSig2CV.loc[~mutSig2CV.index.duplicated(keep='first')]
sigPAMs = list(mutSig2CV.index[mutSig2CV['q']<=0.05])

# Get list of significantly CNA amplified genes
ampGenes = []
ampLoci = {}
amp1 = pd.read_csv(params['gistic_path']+'/'+params['amp_path'],index_col=0,sep='\t')
for col1 in amp1.columns:
    if float(amp1[col1]['residual q value'])<=0.05:
        ampGenes += [i.lstrip('[').rstrip(']') for i in list(amp1[col1].dropna()[3:])]
        ampLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))

# Get list of significantly CNA deleted genes
delGenes = []
delLoci = {}
del1 = pd.read_csv(params['gistic_path']+'/'+params['del_path'],index_col=0,sep='\t')
for col1 in del1.columns:
    if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
        delGenes += [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])]
        delLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))

# Convert gene ids for amp and del genes
ampGenes = [n1.loc[i] for i in ampGenes if i in n1.index]
delGenes = [n1.loc[i] for i in delGenes if i in n1.index]

# Load up somatically mutated genes
somMuts = pd.read_csv(params['som_mut_path'],index_col=0,header=0)
if not somMuts.index.dtype=='int64':
    somMuts = somMuts.loc[somMuts.index.map(lambda x: x in n1.index)]
    somMuts.index = somMuts.index.map(lambda x: n1.loc[x])
somMuts = somMuts.loc[~somMuts.index.duplicated(keep='first')]

# Read in gistic all_data_by_genes file
with open(params['gistic_path']+'/'+params['thresh_path'],'r') as inFile:
    tmp = inFile.readline().strip().split('\t')
    numCols1 = len(tmp)
d1 = pd.read_csv(params['gistic_path']+'/'+params['thresh_path'],index_col=0,sep='\t').drop(tmp[1], axis = 1)
d1.columns = [i[:12] for i in d1.columns]
d1 = d1.loc[d1.index.map(lambda x: x.split('|')[0] in n1.index)]
d1.index = d1.index.map(lambda x: n1.loc[x.split('|')[0]])
d1.index.name = 'Locus ID'

# Removing sex chromosomes (issues in CNA analysis) from d1
lociThresh = d1['Cytoband']
include = []
for i in lociThresh:
    if not(i[0]=='X' or i[0]=='Y'):
        include.append(True)
    else:
        include.append(False)
d1 = d1.loc[lociThresh[include].index].drop('Cytoband', axis = 1)

# Removing sex chromosomes from ampLoci
delMes = [i for i in ampLoci if i[0]=='X' or i[0]=='Y']
for delMe in delMes:
    del ampLoci[delMe]

# Removing sex chromosomes from delLoci
delMes = [i for i in delLoci if i[0]=='X' or i[0]=='Y']
for delMe in delMes:
    del delLoci[delMe]

# Make sure somMuts and gistic have same samples
somMuts = somMuts[list(set(d1.columns).intersection(somMuts.columns))]
d1 = d1[list(set(d1.columns).intersection(somMuts.columns))]

# Get rid of duplicated rows
d1 = d1[~d1.index.duplicated(keep='first')]
min_mut_freq = params['min_mut_freq']
perm_pv = params['perm_pv']

# Print out some useful information
print('\tSize of CNA matrix: '+str(d1.shape))
print('\tSize of somatic mutation matrix: '+str(somMuts.shape))
print('Finished loading data.')

# Make the output directory if it doesn't exists already
if not os.path.exists(params['output_path']):
    os.mkdir(params['output_path'])


########################
## Begin onocoMerging ##
########################
# Cutoff somatic mutations based on the minimum mutation frequency (mf)
freq1 = somMuts.sum(axis=1)/len(list(somMuts.columns))
somMutPoint = freq1[freq1>=params['min_mut_freq']].index

# Precompute positive and negative dichotomized matrices
posdicot = (lambda x: 1 if x>=2 else 0)
posD1 = d1.applymap(posdicot)
negdicot = (lambda x: 1 if x<=(-2) else 0)
negD1 = d1.applymap(negdicot)
print('Finished precomputing dichotomized matrices.')

# Merge loci for deletions and amplifications
lociCNAgenes = {}
lociCNA = pd.DataFrame(columns=d1.columns)
print('Combining deletion loci...')
for loci1 in delLoci:
    # Get matrix of CNAs for genes in loci
    dt = negD1.loc[delLoci[loci1]]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAdel'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]

print('Combining amplification loci..')
for loci1 in ampLoci:
    # Get matrix of CNAs for genes in loci
    dt = posD1.loc[ampLoci[loci1]]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAamp'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]


# Calculating Shallow Coincidence
calcShallowCoincidence = ()
shallowCoincidence = {'Act':{},'LoF':{}}
totalPats = len(d1.columns)

# Activating Shallow Coincidence
for loci1 in ampLoci.keys():
    for ampGene in ampLoci[loci1]:
        if ampGene in somMuts.index and ampGene in d1.index:
            ampProf = d1.loc[ampGene]
            ampPats = set(ampProf.loc[ampProf.apply(np.sign) == 1].index)
            somProf = somMuts.loc[ampGene]
            somPats = somProf.loc[somProf == 1].index
            coincPats = ampPats.intersection(somPats)
            shallowCoincidence['Act'][ampGene] = len(coincPats)/totalPats

# Loss of Function Shallow Coincidence
for loci1 in delLoci.keys():
    for delGene in delLoci[loci1]:
        if delGene in somMuts.index and delGene in d1.index:
            delProf = d1.loc[delGene]
            delPats = set(delProf.loc[delProf.apply(np.sign) == -1].index)
            somProf = somMuts.loc[delGene]
            somPats = somProf.loc[somProf == 1].index
            coincPats = delPats.intersection(somPats)
            shallowCoincidence['LoF'][delGene] = len(coincPats)/totalPats

# Make combined matrix
# LoF = deletions + somatic point mutations
# Act = amplifications + somatic point mutations
print('Starting somatic mutations...')
pamLofAct = {}
freq = {}
for s1 in somMutPoint:
    if s1>0:
        if not str(s1) in pamLofAct:
            pamLofAct[str(s1)] = {}
        tmpSom = somMuts.loc[s1]
        # If potential PAM, store PAM
        if not (str(s1)+'_PAM' in pamLofAct[str(s1)] or sum(tmpSom)==0):
            pamLofAct[str(s1)][str(s1)+'_PAM'] = tmpSom
        if (s1 in negD1.index and s1 in posD1.index):
            tmpNeg = negD1.loc[s1]
            tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
            tmpLoF[tmpLoF==2] = 1
            tmpPos = posD1.loc[s1]
            tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
            tmpAct[tmpAct==2] = 1
            if not s1 in freq:
                freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
        else:
            if not s1 in freq:
                freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}

print('Starting deletions...')
for loci1 in delLoci:
    for s1 in set(delLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Lof
            if s1 in somMuts.index and s1 in negD1.index:
                tmpSom = somMuts.loc[s1]
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                tmpLoF[tmpLoF==2] = 1
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                tmpAct[tmpAct==2] = 1
                if not s1 in freq:
                    freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
                # Store LoF
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_LoF' in pamLofAct[str(s1)] or tmpLoF.equals(tmpSom) or tmpLoF.equals(tmpNeg)):
                    pamLofAct[str(s1)][str(s1)+'_LoF'] = tmpLoF

print('Starting amplifications...')
for loci1 in ampLoci:
    for s1 in set(ampLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Act
            if s1 in somMuts.index and s1 in posD1.index:
                tmpSom = somMuts.loc[s1]
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                tmpLoF[tmpLoF==2] = 1
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                tmpAct[tmpAct==2] = 1
                if not s1 in freq:
                    freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
                # Store Act
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_Act' in pamLofAct[str(s1)] or tmpAct.equals(tmpSom) or tmpAct.equals(tmpPos)):
                    pamLofAct[str(s1)][str(s1)+'_Act'] = tmpAct

print('Screening for frequency...')
keepPAM = []
keepers = {}
calcSig = []
for s1 in pamLofAct:
    if s1 in freq:
        freqPAM = freq[s1]['PAM']
        freqNeg = freq[s1]['CNAdel']
        freqLoF = freq[s1]['LoF']
        freqPos = freq[s1]['CNAamp']
        freqAct = freq[s1]['Act']
        if freqLoF>=0.05 or freqAct>=0.05 or freqPAM>=0.05:
            print('\t'+''.join([str(i) for i in [n1.index[n1==int(s1)][0]+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
        if freqPAM>0 and freqPAM>=params['min_mut_freq'] and int(s1) in somMutPoint and int(s1) in sigPAMs:
            keepers[str(s1)+'_PAM'] = pamLofAct[str(s1)][str(s1)+'_PAM']
            keepPAM.append(str(s1)+'_PAM')
        if str(s1)+'_LoF' in pamLofAct[str(s1)]  and freqLoF>freqPAM and freqLoF>=params['min_mut_freq'] and shallowCoincidence['LoF'][int(s1)] >= params['min_coinc']:
            keepers[str(s1)+'_LoF'] = pamLofAct[str(s1)][str(s1)+'_LoF']
            calcSig.append(str(s1)+'_LoF')
        if str(s1)+'_Act' in pamLofAct[str(s1)] and freqAct>freqPAM and freqAct>=params['min_mut_freq'] and shallowCoincidence['Act'][int(s1)] >= params['min_coinc']:
            keepers[str(s1)+'_Act'] = pamLofAct[str(s1)][str(s1)+'_Act']
            calcSig.append(str(s1)+'_Act')

## Write out LoF and Act significance file
print('Permutation anlaysis...')
numPermutes = 1000

# Permute to get frequency
def singlePermute(somMutsMF, somCNAsMF):
    tmp1 = pd.Series(np.random.permutation(somMutsMF), index=somMutsMF.index)
    tmp2 = pd.Series(np.random.permutation(somCNAsMF), index=somCNAsMF.index)
    subset1 = set(somMutsMF.index).intersection(somCNAsMF.index)
    return list(tmp1.loc[subset1]+tmp2.loc[subset1])

# Deletions
print('\tPermuting deletions...')
permMF_neg = []
somMutsMF = somMuts.transpose().mean()
somCNAsMF = negD1.transpose().mean()
for i in range(numPermutes):
    permMF_neg += singlePermute(somMutsMF, somCNAsMF)

# Amplifications
print('\tPermuting amplifications...')
permMF_pos = []
somMutsMF = somMuts.transpose().mean()
somCNAsMF = posD1.transpose().mean()
for i in range(numPermutes):
    permMF_pos += singlePermute(somMutsMF, somCNAsMF)

# Write Permutation Analysis file Lof_Act_sig
lofActSig = pd.DataFrame(columns = ['Symbol', 'Type','Freq','Emp.p_value'], index = calcSig)
for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_LoF'))][0]
        lofActSig['Type'].loc[sig1] = 'LoF'
        lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_LoF')]['LoF']
    elif sig1.find('Act')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_Act'))][0]
        lofActSig['Type'].loc[sig1] = 'Act'
        lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_Act')]['Act']

permdict1 = {}
permdict1['LoF'] = {}
permdict1['Act'] = {}
oMfreqs = [freq[sig1[:-4]][sig1.split('_')[1]] for sig1 in calcSig]
for f in set(oMfreqs):
    permdict1['LoF'][f] = float(len([i for i in permMF_neg if i >= f]))/len(permMF_neg)
    permdict1['Act'][f] = float(len([i for i in permMF_pos if i >= f]))/len(permMF_pos)

for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['LoF'][freq[sig1.rstrip('_LoF')]['LoF']]
    elif sig1.find('Act')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['Act'][freq[sig1.rstrip('_Act')]['Act']]

if len(lofActSig)>0:
    lofActSig['q_value'] = multipletests(lofActSig['Emp.p_value'], 0.05, method='fdr_bh')[1]
    lofActSig.sort_values('q_value').to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    # Screen out LoF and Act that don't meet significance cutoffs
    keepLofAct = list(lofActSig.index[lofActSig['q_value']<=params['perm_pv']])
else:
    lofActSig.to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    keepLofAct = []

# Screen out PAMs that are LoF/Act
newKeepPAM = []
for pam1 in keepPAM:
    found = 0
    tmp1 = pam1.split('_')[0]
    for lofAct in keepLofAct:
        if tmp1==lofAct.split('_')[0]:
            found = 1
    if found==0:
        newKeepPAM.append(pam1)

## Screen out loci that have a representative gene
# Mutations that are at or above minimum mutation frequency cutoff
highFreqLoci = lociCNA[lociCNA.mean(axis=1)>=params['min_mut_freq']]

# Figure out what loci are explained by current Act or LoF genes
explainedLoc = []
for locus1 in highFreqLoci.index:
    genesInLocus = [i for i in lociCNAgenes[locus1] if (str(i)+'_Act' in keepLofAct or str(i)+'_LoF' in keepLofAct)]
    if len(genesInLocus)>0:
        explainedLoc.append(locus1.split('_')[0])

# Screen out all other loci in that region
keepLoc = []
for locus1 in highFreqLoci.index:
    if not locus1.split('_')[0] in explainedLoc:
        keepLoc.append(locus1)

keepLoc_dict = {}
for locus1 in set(['_'.join([i.split('_')[0],i.split('_')[-1]]) for i in keepLoc]):
    locus2 = locus1.split('_')[0]
    if locus2 in ampLoci.keys():
        keepLoc_dict[locus1] = posD1.loc[ampLoci[locus2]]
    if locus2 in delLoci.keys():
        keepLoc_dict[locus1] = negD1.loc[delLoci[locus2]]

def mode2(pat1col):
    tmpser = pat1col.mode()
    if len(tmpser)>1:
        tmp10 = 1
    else:
        tmp10 = tmpser.iloc[0]
    return tmp10

keepLoc_df = pd.DataFrame(columns = d1.columns)
tmpMode = pd.Series(index = d1.columns)
for locus1 in keepLoc_dict.keys():
    for pat1 in keepLoc_dict[locus1].columns:
         tmpMode.loc[pat1] = mode2(keepLoc_dict[locus1][pat1])
    if tmpMode.mean()>=0.05:
        keepLoc_df.loc[locus1] = tmpMode

keepLoc_df = keepLoc_df.applymap(int)

####################################
## Compile OncoMerge output files ##
####################################
finalMutFile = pd.concat([pd.DataFrame(keepers).transpose().loc[newKeepPAM].sort_index(), pd.DataFrame(keepers).transpose().loc[keepLofAct].sort_index(), keepLoc_df.sort_index()], sort=True)

# Rename all loci with only one gene
ind_list = finalMutFile.index.tolist()
for locus1 in keepLoc_df.index:
    splitUp = locus1.split('_')
    if splitUp[-1]=='CNAamp' and len(ampLoci[splitUp[0]])==1:
        idx = ind_list.index(locus1)
        ind_list[idx] = str(ampLoci[splitUp[0]][0]) + '_CNAamp'
    if splitUp[-1]=='CNAdel' and len(delLoci[splitUp[0]])==1:
        idx = ind_list.index(locus1)
        ind_list[idx] = str(delLoci[splitUp[0]][0]) + '_CNAdel'

finalMutFile.index = ind_list

# Write file
finalMutFile.to_csv(params['output_path']+'/oncoMerge_mergedMuts.csv')

## Write out loci
# Prepare for writing out
writeLoci = ['Locus_name,Genes']
for locus1 in keepLoc_df.index:
    splitUp = locus1.split('_')
    if splitUp[-1]=='CNAamp':
        writeLoci.append(locus1+','+' '.join([str(i) for i in ampLoci[splitUp[0]]]))
    if splitUp[-1]=='CNAdel':
        writeLoci.append(locus1+','+' '.join([str(i) for i in delLoci[splitUp[0]]]))

# Write out file
with open(params['output_path']+'/oncoMerge_CNA_loci.csv','w') as outFile:
    outFile.write('\n'.join(writeLoci))

print('Done.')

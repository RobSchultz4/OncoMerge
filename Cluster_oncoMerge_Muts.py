# Correlation Analysis
# Correlate the mutational profiles of Muts within each Cancer. Print out 32 heat maps and clustered oncoMerge outputs



import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as spc
import numpy as np
import os
import pickle as pkl
import argparse
parser = argparse.ArgumentParser(description='Used to choose a cancer type')
parser.add_argument('-o', '--output_path', help='Path to output Location', type = str)
args = parser.parse_args()

def FindLoci(inp,cancer):
    gene1, mutType = inp.split('_')
    if 'p' in gene1 or 'q' in gene1:
        ret = gene1
        LociSet = {}
    elif mutType == 'Act' or mutType == 'CNAamp':
        LociSet = ins[cancer]['ampLoci']
    elif mutType == 'LoF' or mutType == 'CNAdel':
        LociSet = ins[cancer]['delLoci']
    else:
        ret = 'NA'
        LociSet = {}
    for loci1 in LociSet.keys():
        if int(gene1) in LociSet[loci1]:
            ret = loci1
    return ret


def ent2sym2(ents):
    if len(ents)>0:
        syms = []
        if isinstance(list(ents)[0], tuple):
            syms = [ent2sym_dict[ent[1]] for ent in ents if ent[1] in ent2sym_dict.keys()]
        if isinstance(list(ents)[0], str):
            syms = [ent2sym_dict[ent] if ent in ent2sym_dict.keys() else ent for ent in ents ]
        if isinstance(list(ents)[0], int):
            syms = [ent2sym_dict[str(ent)] for ent in OMgenes[cancer] if str(ent) in ent2sym_dict.keys()]
        syms = set(syms)
    else:
        syms = list()
    return ' | '.join(syms)
#outs = pd.read_pickle('outs.pkl')


print('PKL FOUND. loading ins and outs from pkl')
ins = pd.read_pickle('output/ins.pkl')
n2 = ins['ACC']['n1'].reset_index().set_index('Locus ID')
n1 = ins['ACC']['n1']
n2.index = n2.index.map(str)
ent2sym_dict = {i:v[0] for i,v in zip(n2.index.values,n2.values)}


cancers = sorted(['UCS', 'SKCM', 'KICH','ESCA','ACC','DLBC','READ','COAD','GBM','LGG','PCPG','BLCA','UCEC','THCA','CESC','THYM','LIHC','CHOL','HNSC','UVM','PAAD','TGCT','LUSC','MESO','OV','SARC','KIRP','STAD','PRAD','LUAD','BRCA','KIRC'])

# Read in oncoMerge mutFiles
omdf = {}
smCor = {}
#cancer = 'ACC'
for cancer in cancers:
    omdf[cancer] = pd.read_csv(args.output_path + 'oncoMerged_'+cancer+'/oncoMerge_mergedMuts.csv', index_col = 0)
    smCor[cancer] = omdf[cancer].T.corr()

dict1 = {}
dict1['Muts'] = {}
dict1['corrMuts'] = {}
pairs = {}
for cancer in cancers:
    pairs[cancer] = {}
    dict1['Muts'][cancer] = len(smCor[cancer].columns)
    dict1['corrMuts'][cancer] = 0
    for mut1 in smCor[cancer].columns:
        pairs[cancer][mut1] = []
        tmpcol = (smCor[cancer][mut1] >= 0.8) & (smCor[cancer][mut1] < 1)
        for i in smCor[cancer].index:
            if tmpcol.loc[i]:
                pairs[cancer][mut1].append(i)

for cancer in cancers:
    for mut1 in pairs[cancer].keys():
        if not pairs[cancer][mut1] == []:
            dict1['corrMuts'][cancer] +=1

dict1['Percent Correlated'] = {}
for cancer in cancers:
    dict1['Percent Correlated'][cancer] = dict1['corrMuts'][cancer]/dict1['Muts'][cancer]

#initalize data structures
clustList = {}
clustDict= {}
BigDists = {}
count = {}
dc = {}
pdist = {}
for cancer in cancers:
    dc[cancer] = omdf[cancer].T.corr().values
    print(len(dc[cancer]))
    pdist[cancer] = spc.distance.pdist(omdf[cancer], 'correlation')
    BigDists[cancer] = []
    clustList[cancer] = []
    clustDict[cancer] = {}
    count[cancer] = 0


# Build Linkage
linkage = {}
for cancer in cancers:
    linkage[cancer] = spc.linkage(pdist[cancer], method='complete')


# Find dists for each cancer
dists= {}
distInd= {}
for cancer in cancers:
    dists[cancer] = sorted([max(i) for i in spc.dendrogram(linkage[cancer])['dcoord']],reverse = True)
    dists[cancer].append(0)
    distInd[cancer] = 0


def getClusterInfo(i):
    tmp0 = omdf[cancer].T.corr()[idx==i]
    tmp01 = tmp0[tmp0.index]
    if sum(idx==i)>1:
        tmp02 = tmp01.values[np.triu_indices(len(tmp01),k=1)]
    else:
        tmp02 = np.array([1])
    clusterInfo = {'size': sum(idx==i), 'minCorr': np.amin(tmp02),'genes':list(tmp01.index)}
    return clusterInfo


# Find cut points for clusters and capture genes in each highly correlated cluster
oops = []
genesAccountedFor ={}
for cancer in cancers:
    print('#####' + cancer + '#####')
    genesAccountedFor[cancer] = []
    clusterInfo = {'NoClusts':{'minCorr':0, 'size':0}}
    while any([clusterInfo[clust1]['minCorr'] < 0.8 for clust1 in clusterInfo.keys()]) and not len(dists[cancer])==distInd[cancer]:
        dist = dists[cancer][distInd[cancer]]
        idx = spc.fcluster(linkage[cancer], dist, 'distance')
        clusterInfo = {}
        # Calculate average correlation within clusters
        for i in set(idx):
            clusterInfo[i] = getClusterInfo(i)
            if clusterInfo[i]['minCorr']>=0.8:
                if all(not gene1 in genesAccountedFor[cancer] for gene1 in clusterInfo[i]['genes']):
                    genesAccountedFor[cancer].extend(clusterInfo[i]['genes'])
                    mutTypes = set([i.split('_')[1] for i in clusterInfo[i]['genes']])
                    #Only needed if we find clusters with this happening
                    if ('LoF' in mutTypes or 'CNAdel' in mutTypes) and ('Act' in mutTypes or 'CNAamp' in mutTypes):
                        oops.append((cancer,count[cancer],clusterInfo[i],i,dist))
                    #then, split into separate clusters
                    clustDict[cancer][count[cancer]] = {}
                    clustDict[cancer][count[cancer]]['minCorr'] = clusterInfo[i]['minCorr']
                    clustDict[cancer][count[cancer]]['genes'] = clusterInfo[i]['genes']
                    clustDict[cancer][count[cancer]]['symbols'] = ent2sym2([i.split('_')[0] for i in clusterInfo[i]['genes']])
                    clustDict[cancer][count[cancer]]['loci'] = [FindLoci(i,cancer) for i in clusterInfo[i]['genes']]
                    clustDict[cancer][count[cancer]]['distance'] = dist
                    clustList[cancer].append(clusterInfo[i]['genes'])
                    count[cancer] += 1
        distInd[cancer] += 1
        len(clustList[cancer])
    #clustList[cancer] = clusterInfo



'''for cancer in cancers:
    for count1 in range(count[cancer]):
        gl1 = clustList[cancer][count1]
        clustDict[cancer][count1] = {}
        clustDict[cancer][count1]['genes'] = gl1
        clustDict[cancer][count1]['symbols'] = ent2sym2(gl1)
        clustDict[cancer][count1]['loci'] = [FindLoci(i, cancer) for i in gl1]
        clustDict[cancer][count1]['distances'] = BigDists[cancer][count1]
'''

# Write out Cluster Membership
with open(args.output_path + 'analysis/SYGNAL_Input_Cluster_Membership.csv', 'w') as outFile:
    outFile.write('Tumor Type, Cluster Name, Cluster Members, Symbols, Loci, Distance, Minimum Correlation\n') #outFile header
    for cancer in cancers:
        for clust1 in clustDict[cancer]:
            line = [cancer,str(clust1),' | '.join(clustDict[cancer][clust1]['genes']), ''.join(clustDict[cancer][clust1]['symbols']),' | '.join(clustDict[cancer][clust1]['loci']),str(clustDict[cancer][clust1]['distance']),str(clustDict[cancer][clust1]['minCorr']) ]
            outFile.write(','.join(line)+'\n')

# Save out the contents of each cluster
with open(args.output_path + 'analysis/SYGNAL_Input_Clustering_Comparison.csv', 'w') as outFile:
    for cancer in cancers:
        ents = [cancer,'Entrez']
        syms = [cancer,'Symbols']
        loci = [cancer,'Loci']
        distances = [cancer,'Distance']
        Correlation= [cancer,'Correlation']
        for clust1 in clustDict[cancer]:
            ents.append(' | '.join(clustDict[cancer][clust1]['genes']))
            syms.append(''.join(clustDict[cancer][clust1]['symbols']))
            loci.append(' | '.join(clustDict[cancer][clust1]['loci']))
            distances.append(str(clustDict[cancer][clust1]['distance']))
            Correlation.append(str(clustDict[cancer][clust1]['minCorr']))
        outFile.write(','.join(ents)+'\n')
        outFile.write(','.join(syms)+'\n')
        outFile.write(','.join(loci)+'\n')
        outFile.write(','.join(distances)+'\n')
        outFile.write(','.join(Correlation)+'\n')
        outFile.write('\n')


# Calculate distance statistics
from statistics import stdev
BigDists = [i for j in [[clustDict[cancer][i]['distance'] for i in clustDict[cancer].keys() if not clustDict[cancer][i]['minCorr'] == 1 ] for cancer in cancers] for i in j]
avgDist = sum(BigDists)/len(BigDists)
stdDist = stdev(BigDists)
minDist = min(BigDists)
maxDist = max(BigDists)

#plt.hist(BigDists)

#plt.show()

def mode2(pat1col):
    tmpser = pat1col.mode()
    if len(tmpser)>1:
        tmp10 = 1
    else:
        tmp10 = tmpser.iloc[0]
    return tmp10

# combine Clusters by taking the frequency of the best feature.
#If tied, choose the ones with the most other genes, otherwise, idk...
mishap = {}
count1 = 0
omdf1 = {}
clusterNames = {}
for cancer in cancers:
    clusterNames[cancer] = {}
    line = []
    omdf1[cancer] = omdf[cancer]
    tmpProf = pd.Series(index = omdf[cancer].columns)
    for clust1 in clustDict[cancer].keys():
        tmpgenes = clustDict[cancer][clust1]['genes']
        tmpFreqs = omdf[cancer].loc[tmpgenes].mean(axis = 1)
        tmpProfInd = tmpFreqs.loc[tmpFreqs == tmpFreqs.max()].index[0]
        tmpProf = omdf[cancer].loc[tmpProfInd]
        if len(tmpgenes) == 1:
            tmpProf.name = tmpgenes[0]
            clusterNames[cancer][clust1] = tmpgenes[0]
        else:
            clusterNames[cancer][clust1] = 'cluster_' + str(clust1) + '_' + tmpProfInd.split('_')[1]
            tmpProf.name = 'cluster_' + str(clust1) + '_' + tmpProfInd.split('_')[1]
        omdf1[cancer] = omdf1[cancer].drop(tmpgenes)
        omdf1[cancer] = omdf1[cancer].append(tmpProf)
        if tmpProf.mean() < 0.05:
            print(str(count1) +'. Yikes: '+cancer+': ############################################')
            print('Mutational Frequencies:')
            omdf[cancer].loc[tmpgenes].mean(axis = 1)
            mishap[count1] = omdf[cancer].loc[tmpgenes]
            count1 += 1

for cancer in cancers:
    omdf1[cancer].to_csv(args.output_path + 'analysis/Clustered_Profiles_minCorr_maxFreq/'+cancer+'oncoMerge_mergedMuts.csv')


dict1['Features'] = {}
dict1['Clusters'] = {}
dict1['Solo Mutations'] = {}
for cancer in cancers:
    dict1['Features'][cancer] = len(clustDict[cancer])
    clusters = 0
    solos = 0
    for clust1 in clustDict[cancer].keys():
        if len(clustDict[cancer][clust1]['genes']) > 1:
            clusters+=1
        else:
            solos += 1
    dict1['Clusters'][cancer] = clusters
    dict1['Solo Mutations'][cancer] = solos

with open(args.output_path + 'analysis/Mut_Correlation_Table.csv','w') as outFile:
    outFile.write('Tumor Type, Mutations,Correlated Mutations,Percent Correlated,Features,Clusters,Solo Mutations\n') # header
    for cancer in cancers:
        line = []
        line.append(cancer)
        line.append(str(dict1['Muts'][cancer]))
        line.append(str(dict1['corrMuts'][cancer]))
        line.append(str(round(dict1['Percent Correlated'][cancer]*100,1)))
        line.append(str(dict1['Features'][cancer]))
        line.append(str(dict1['Clusters'][cancer]))
        line.append(str(dict1['Solo Mutations'][cancer]))
        outFile.write(','.join(line)+'\n')


#freq_vote_1_3 = [len([i for i in mishap[m].sum() if i>=len(mishap[m])/3])/len(mishap[m].columns) for m in mishap.keys()]


# CESC, CESC, CESC, CESC, CESC, HNSC, HNSC, LGG, LIHC
# 80.2, 82.8, 81.4, 80.3, 83.6

'''
spc.dendrogram(linkage)
plt.show()
'''











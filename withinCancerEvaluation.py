# oncoMerge Cancer Specific Plots
import json
import argparse
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.')
parser.add_argument('-cf', '--config_path', help='Path to JSON encoded configuration file, overrides command line parameters', type = str)
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default = '.')
args = parser.parse_args()


params = args.__dict__
if args.config_path:
    with open(args.config_path, "r") as cfg:
        tmp = json.loads(cfg.read())
        for i in tmp:
            params[i] = tmp[i]


summaryMatrix = pd.read_csv(params['output_path']+'/oncoMerge_summaryMatrix.csv')

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

#Plots Violins showing distribution of mutation frequency for different Mutation Types.
with PdfPages(args.output_path+ '/FreqViolins_byMutType_swarm.pdf') as pdf:
    data = summaryMatrix['Final_freq','Final_mutation_type']
    ax = sns.violinplot(x = 'Final_mutation_type', y = 'Final_freq', scale = 'count', inner = None, order = ['PAM','Act','CNAamp','LoF','CNAdel'], data = data)
    ax = sns.swarmplot( x = 'Final_mutation_type', y = 'Final_freq', color = 'k',     size = 2,     order = ['PAM','Act','CNAamp','LoF','CNAdel'], data = data) #, ax = ax)
    plt.ylim(float(mf), 1)
    plt.title('Frequency of Mutation'
    pdf.savefig()
    plt.close()


# Create a table showing the Number of recovered Genes
recoveredMuts_df = summaryMatrix.loc[(summaryMatrix['PAM_freq'] < 0.05) & (summaryMatrix['CNA_freq'] < 0.05)]
recovered_MutTypes_Raw_df = pd.DataFrame(index = ['Recovered_mutations'], columns = ['Act','LoF'])
for sm in ['Act', 'LoF']:
        recovered_MutTypes_Raw_df.loc['Recovered_mutations',sm] = sum(recoveredMuts_df['Final_mutation_type'] == sm)

recovered_MutTypes_Raw_df.to_csv(args.output_path + '/recovered_MutTypes_Raw.csv')

# Create a Stacked Bar Plot showing the extension of the mutational frequency by including CNA data.
with PdfPages(args.output_path+'/StackedBars.pdf') as pdf:
    summaryMatrix = summaryMatrix.loc[list(summaryMatrix['Final_mutation_type'].sort_values().index)]
    N = len(summaryMatrix)
    ind = np.arange(N)
    p1 = plt.bar(ind,summaryMatrix['PAM_only'], color = 'y')
    p2 = plt.bar(ind,summaryMatrix['Coincidence_rate'], bottom = summaryMatrix['PAM_only'] , color = 'g')
    p3 = plt.bar(ind,summaryMatrix['CNA_only'], bottom =summaryMatrix['Coincidence_rate']+ summaryMatrix['PAM_only'], color = 'b')
    p4 = plt.plot([-1,len(ind)], [0.05,0.05], 'r', label='Mutation Frequency Threshold', linewidth = 1/len(ind))
    plt.ylabel('Frequency')
    plt.title('oncoMerged Mutation Frequencies')
    plt.xticks(ind,tuple(summaryMatrix['Symbol']), rotation = 90)
    plt.legend(tuple([p1[0],p2[0],p3[0],p4[0]]), tuple(['PAM only', 'Both','CNA only', 'Threshold']))
    plt.tight_layout()
    pdf.savefig()
    plt.close()







###################
#Import packages
###################
import pandas as pd
import pickle
import random
import numpy as np

import sys
print(sys.executable)
print(sys.version)
print(sys.path)
sys.path.extend(['/scratch/research/repos', 
                 '/scratch/research/repos/ness_fasta', 
                 '/scratch/research/repos/annotation', 
                 '/scratch/research/repos/ness_fastq', 
                 '/scratch/research/repos/ness_vcf', 
                 '/scratch/research/repos/vcf2fasta'])
from ness_vcf import SFS

##################
#Define functions
##################

def bootstrap_with_replacement(sample_list, size, rep):

    '''randomly sample n elements from sample_list using rep to set random.seed'''
    
    random.seed(rep)
    sample = np.random.choice(sample_list,size)
    return sample

def get_total_SFSs (list_of_indexes, SFS_dataframe):

    '''calculates the total SFS vector of a subset list of indexes \
    given a dataframe that contains selected_SFS and neutral_SFS'''
    
    selected = [sum(i) for i in zip(*SFS_dataframe.loc[list_of_indexes, "selected_SFS"])]
    neutral = [sum(i) for i in zip(*SFS_dataframe.loc[list_of_indexes, "neutral_SFS"])]
    
    return selected, neutral

####################
#Define variables
####################

n_sample = int(snakemake.params.sample)
#pool = snakemake.input.pool
df = pd.read_pickle(str(snakemake.input))
sample_list = df.index

####################
#Random sampling
####################

sample = bootstrap_with_replacement(sample_list = sample_list,
                                    size = snakemake.config['RSAMPLE_SIZE'],
                                    rep = n_sample)

selected, neutral = get_total_SFSs(list_of_indexes = sample, 
                                   SFS_dataframe = df)

selected = SFS(selected).fold().sfs
neutral = SFS(neutral).fold().sfs

selected_div = [str(i) for i in [1,sum(df.sites0),sum(df.diffs0)]]
neutral_div = [str(i) for i in [0, sum(df.sites4), sum(df.diffs4)]]

dictionary = {'selected': selected,
             'neutral': neutral,
             'selected_div': selected_div,
             'neutral_div': neutral_div}
             
with open(str(snakemake.params.prefix), "wb") as f:
    pickle.dump(dictionary, f)
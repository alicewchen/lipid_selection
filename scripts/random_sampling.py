###################
#Import packages
###################
import pandas as pd
import pickle
import random
import numpy as np
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
dictionary = {'selected': selected,
             'neutral': neutral}
             
with open(str(snakemake.params.prefix), "wb") as f:
    pickle.dump(dictionary, f)
#####################################
# Import packages and modules
#####################################
import pickle
import sys
sys.path.extend(['/scratch/research/repos', 
                 '/scratch/research/repos/ness_fasta', 
                 '/scratch/research/repos/annotation', 
                 '/scratch/research/repos/ness_fastq', 
                 '/scratch/research/repos/ness_vcf', 
                 '/scratch/research/repos/vcf2fasta'])
from ness_vcf import SFS

#####################
# Define variables
#####################

dictionary = pickle.load(open(str(snakemake.input), "rb"))
sfs_input_path = str(snakemake.output[0])
divergence_file_path = str(snakemake.output[1])


##########################
# Generate SFS_input.txt
##########################

neutral = [str(i) for i in SFS(dictionary['neutral']).fold().sfs]
selected = [str(i) for i in SFS(dictionary['selected']).fold().sfs]

neutral = " ".join(neutral)
selected = " ".join(selected)

M = str(1)
MIN_ALLELE = str(18)
ls = [M, MIN_ALLELE, selected, neutral]

with open(sfs_input_path, 'w') as f:
    f.write("\n".join(ls))


###############################     
# Generate divergence_file.txt
############################### 

with open(divergence_file_path, 'w') as f:
    f.write(" ".join(dictionary['selected_div'])+"\n") 
    f.write(" ".join(dictionary['neutral_div']))
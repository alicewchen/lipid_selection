import pandas as pd
import pickle as pk

with open(str(snakemake.output), "w") as f: 
    with open(str(snakemake.input[0]), "rb") as pool_file:
        df = pd.read_pickle(pool_file)
        output = [str(snakemake.input[0])+"\t"+str(len(df))+"\n"]
    with open(str(snakemake.input[1]), "rb") as pool_file:
        df = pd.read_pickle(pool_file)
        output.append(str(snakemake.input[1])+"\t" + str(len(df)))
    f.writelines(output)
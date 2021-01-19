with open(str(snakemake.input), "r") as f:
    input = f.read()
with open("../data/workflow/bootstrap/"+str(snakemake.params.pool)+"_bootstrap_values.txt", "a") as f:
    f.write(input)
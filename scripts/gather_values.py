ls = snakemake.input
for i in ls:
    with open(str(i), "r") as f:
        input = f.read()
    with open(str(snakemake.output), "a") as f:
        f.write(input)
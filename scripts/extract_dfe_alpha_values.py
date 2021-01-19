import re
with open(str(snakemake.input.omega_output)) as f:
    file = f.read()
    match = re.search(r"Fixation prob of deleterious mutation (\d*\.?\d*)", file)
    pFDM = match.group(1)
    match = re.search(r"alpha (\d*\.?\d*) omega_A (\d*\.?\d*)" ,file)
    alpha = match.group(1)
    omega = match.group(2)
with open(str(snakemake.input.prop_output)) as f:
    file = f.read()
    match = re.search(r"0.000000 1.000000 (\d*\.?\d*) 1.000000 10.000000 (\d*\.?\d*) 10.000000 100.000000 (\d*\.?\d*) 100.000000 -99.000000 (\d*\.?\d*) ", file)
    NeS_0_1, NeS_1_10, NeS_10_100, NeS_100_99 = match.group(1), match.group(2), match.group(3), match.group(4)
    
ls = [alpha, omega, pFDM, NeS_0_1, NeS_1_10, NeS_10_100, NeS_100_99]
with open(str(snakemake.output), 'w') as f:
    f.write("\t".join(ls)+"\n")
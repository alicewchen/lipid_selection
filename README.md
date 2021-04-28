# Lipid selection

## Description
This project estimates selection of candidate genes responsible for lipid content in *C. reinhardtii*. 

## Getting started

Set project folder `lipid_selection/` as the working directory.

### Install snakemake environment
```
conda env create -f env/snakemake.yml
```

To check the snakemake installation, run the following:
```
cd test_workflow
snakemake --use-conda --cores 1
```

### Install the conda envrionment used to run the rest of the python codes
```
conda env create -f env/main.yml
```

## Running the DFE-𝛼 snakemake workflow

Set working directory as `lipid_selection/workflow`.
```
cd workflow
```

### Perform dry run of Snakemake
```
snakemake -n
```

### Run snakemake
```
snakemake --use-conda --cores 10
```
`--cores` option selects the number of cores used.
If snakemake cannot find DFE-alpha commands, run something like the following to locate the DFE-alpha source folder.
```
export PATH=$PATH:/scratch/research/tmp_apps/dfe-alpha-release-2.16/
source ~/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
```


### Re-run with parameter changes in the config file
Sometimes changing a config file affects only a portion of the workflow, so you don't want to re-run the entire workflow again. Running the following command only changes the output file affected by the config change. If an output file has already been generated previously, it will not be updated to reflect the config change. 
```
snakemake -R snakemake --list-params-changes --cores 10
```

### Remove all snakemake outputs
Snakemake doesn't rerun if the output files are already generated, so remove all the snakemake outputs first. (You don't need a lot of cores to do this quickly)
```
snakemake --delete-all-output --cores 2
```


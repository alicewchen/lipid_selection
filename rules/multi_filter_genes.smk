subworkflow multi_filter_genes:
    workdir:
        "../path/to/otherworkflow"
    snakefile:
        "../path/to/otherworkflow/Snakefile"
    configfile:
        "path/to/custom_configfile.yaml"

rule a:
    input:
        otherworkflow("test.txt")
    output: ...
    shell:  ...
    
rule filter_genes:
    input:
        "../data/workflow/clean_primary_data.pk",
        "../data/workflow/sampled_genes.pk",
        "../data/workflow/SFS_and_divergence.pk"
    output:
        temp("../data/workflow/pools/{pool}.pk")
    params:
        prefix = "../data/workflow/pools/{pool}.pk"
    config: 
        "../config/candidate_config.yaml"
    log:
        notebook="log/notebook/01_Filter_genes_processed_{pool}.ipynb"
    notebook:
        "../analysis/04_filter_data/01_Filter_genes.ipynb"
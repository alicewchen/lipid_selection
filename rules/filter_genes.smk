rule filter_genes:
    input:
        "../data/workflow/clean_primary_data.pk",
        "../data/workflow/sampled_genes.pk",
        "../data/workflow/SFS_and_divergence.pk"
    output:
        temp("../data/workflow/pools/{pool}.pk")
    params:
        prefix = "../data/workflow/pools/{pool}.pk",
    log:
        notebook=config["SM_OUTPUT_DIR"]+"/log/notebook/01_multi_Filter_genes_processed_{pool}.ipynb"
    notebook:
        "../analysis/04_filter_data/01_multi_Filter_genes.ipynb"
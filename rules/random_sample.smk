rule random_sample:
    input:
        rules.filter_genes.output
    output: 
        temp("../data/workflow/samples/{pool}_{sample}.pk")
    params:
        sample = "{sample}",
        prefix = "../data/workflow/samples/{pool}_{sample}.pk",
        candidate_pool = "../data/workflow/pools/candidate.pk"
    script:
        "../scripts/random_sampling.py"
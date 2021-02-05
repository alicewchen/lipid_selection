rule pool_size:
    input:
        expand("../data/workflow/pools/{pool}.pk", pool = pool)
    output:
        expand("{output_dir}/pool_size.txt", output_dir = config["SM_OUTPUT_DIR"])
    params:
        prefix = config["SM_OUTPUT_DIR"]
    script:
        "../scripts/pool_size.py"
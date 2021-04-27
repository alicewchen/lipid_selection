##### Get important values ##### 

rule extract_values:
    input:
        omega_output = rules.run_est_alpha_omega.output,
        prop_output = rules.run_prop_muts_in_s_ranges.output.output_file,
    output:
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/important_values.txt")
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}/",
        sample = "{sample}",
        pool = "{pool}"
    script:
        "../scripts/extract_dfe_alpha_values.py"
        
rule aggregate:
    input:
        expand("../data/workflow/bootstrap/run_{sample}/{pool}/important_values.txt", sample = sample, pool = pool)
    output:
        expand("{prefix}/bootstrap_values.txt", prefix = config["SM_OUTPUT_DIR"])
    params:
        prefix = config["SM_OUTPUT_DIR"]
    script:
        "../scripts/gather_values.py"

#### Create copy of workflow config and output files ####
        
rule save_workflow:
    input:
        "../config/config.yaml",
        "../env/environment.yaml",
        "Snakefile"
    output:
        "{output_dir}/config.yaml",
        "{output_dir}/environment.yaml",
        "{output_dir}/Snakefile"
    params:
        output_dir = config["SM_OUTPUT_DIR"],
    shell:
        "cp {input} {wildcards.output_dir}/"
        
rule export_workflow_DAG:
    output:
        "{output_dir}/snakemake_DAG.pdf"
    params:
        output_dir = config["SM_OUTPUT_DIR"]
    shell:
        "snakemake --dag | dot -Tpdf > {output} "
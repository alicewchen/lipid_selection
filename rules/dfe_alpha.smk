rule generate_dfe_alpha_input:
    input:
        rules.random_sample.output
    output:
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/SFS_input.txt"),
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/divergence_file.txt")
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}"
    group: "dfe_alpha_group"
    script:
        "../scripts/generate_dfe_alpha_input.py"

rule generate_dfe_alpha_config:
    input:
        rules.generate_dfe_alpha_input.output
    output:
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/est_dfe_neut_config.txt"),
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/est_dfe_sel_config.txt"),
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/est_alpha_omega_config.txt")
    params:
        config_template_dir = "../data/workflow/dfe_alpha_config_template/",
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}",
        sample = "{sample}",
        pool = "{pool}"
    group: "dfe_alpha_group"
    script:
        "../scripts/generate_dfe_alpha_config.py"
        
rule run_est_dfe_neutral:
    input:
        rules.generate_dfe_alpha_config.output
    output:
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/run_est_dfe_neutral.log")
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}/",
        sample = "{sample}",
        pool = "{pool}"
    group: "dfe_alpha_group"
    shell:
        "est_dfe -c {input[0]} &> {output} "

rule run_est_dfe_selected:
    input:
        rules.generate_dfe_alpha_config.output,
        rules.run_est_dfe_neutral.output
    output:
        "../data/workflow/bootstrap/run_{sample}/{pool}/run_est_dfe_selected.log"
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}/",
        sample = "{sample}",
        pool = "{pool}"
    group: "dfe_alpha_group"
    shell:
        "est_dfe -c {input[1]} &> {output}"
        
rule run_est_alpha_omega:
    input:
        rules.generate_dfe_alpha_config.output,
        rules.run_est_dfe_selected.output
    output:
        temp("../data/workflow/bootstrap/run_{sample}/{pool}/run_est_alpha_omega.log")
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}/",
        sample = "{sample}",
        pool = "{pool}"
    group: "dfe_alpha_group"
    shell:
        "est_alpha_omega -c {input[2]} &> {output}"
        
rule run_prop_muts_in_s_ranges:
    input:
        rules.run_est_alpha_omega.output
    output:
        output_file = temp("../data/workflow/bootstrap/run_{sample}/{pool}/output_file"),
        output_log = temp("../data/workflow/bootstrap/run_{sample}/{pool}/run_prop_muts_in_s_ranges.log")
    params:
        prefix = "../data/workflow/bootstrap/run_{sample}/{pool}/",
        sample = "{sample}",
        pool = "{pool}"
    group: "dfe_alpha_group"
    shell:
        "prop_muts_in_s_ranges -c {params.prefix}selected/est_dfe.out -o {output.output_file} &> {output.output_log}"
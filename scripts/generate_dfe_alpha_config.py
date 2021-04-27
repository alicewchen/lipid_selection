SFS_input_path = snakemake.input[0]
divergence_file_path = snakemake.input[1]

#est_dfe_neut_config.txt
est_dfe_results_neut_dir_path = str(snakemake.params.prefix)+"/neutral/"

#est_dfe_sel_config.txt
est_dfe_results_sel_dir_path = str(snakemake.params.prefix)+"/selected/"
est_dfe_demography_results_file_path = str(snakemake.params.prefix)+"/neutral/est_dfe.out"

#est_alpha_omega_config.txt
est_dfe_results_file_path =  est_dfe_results_sel_dir_path+"est_dfe.out"
neut_egf_file_path = est_dfe_results_neut_dir_path+"neut_egf.out"
sel_egf_file_path = est_dfe_results_sel_dir_path+"sel_egf.out"
est_alpha_omega_results_file_path = str(snakemake.params.prefix)+"/est_alpha_omega.out"

#est_dfe_neut_config.txt

template = str(snakemake.params.config_template_dir)+"est_dfe_neut_config.txt"

with open(template, 'r') as f:
    f = f.read()
    f = f.replace('SFS_input_path',SFS_input_path)
    f = f.replace('est_dfe_results_neut_dir_path', est_dfe_results_neut_dir_path)
    
with open(str(snakemake.output[0]), "w") as f2:
    f2.write(f)

    
#est_dfe_sel_config.txt

template = str(snakemake.params.config_template_dir)+"est_dfe_sel_config.txt"

with open(template, 'r') as f:
    f = f.read()
    f = f.replace('SFS_input_path',SFS_input_path)
    f = f.replace('est_dfe_results_sel_dir_path', est_dfe_results_sel_dir_path)
    f = f.replace('est_dfe_demography_results_file_path', est_dfe_demography_results_file_path)
    
with open(str(snakemake.output[1]), "w") as f2:
    f2.write(f)
    

    
#est_alpha_omega_config.txt

template = str(snakemake.params.config_template_dir)+"est_alpha_omega_config.txt"

with open(template, 'r') as f:
    f = f.read()
    f = f.replace('divergence_file_path',divergence_file_path)
    f = f.replace('est_alpha_omega_results_file_path',est_alpha_omega_results_file_path)
    f = f.replace('est_dfe_results_file_path', est_dfe_results_file_path)
    f = f.replace('neut_egf_file_path',neut_egf_file_path)
    f = f.replace('sel_egf_file_path',sel_egf_file_path)

with open(str(snakemake.output[2]), "w") as f2:
    f2.write(f)
    

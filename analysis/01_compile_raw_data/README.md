# Extract relevant genes and generate files in intermediate_data_01/

## Workflow

Run the notebooks in order of the file names. 
        
## Important notes:

* Omit 'Blaby_2013_DS8.xlsx' because it is RNA expression data not transcript expression data.

## Notes about files in this folder:

1. 20200716_sampled_genes_key.csv contains manually searched transcript id's for gene symbols that are not found in Phytozome12 using 05_Get_v5.5names_for_experimental_studies.ipynb

2. Use 20200716_sampled_genes_key.csv to get the proper transcript id for each gene symbol mentioned in Summary of primary_data

`04_Data_wrangling...ipynb`:
output: `../../data/intermediate_data_02/sampled_genes.pk`
output columns: 'num_detected', 'num_manipulated', 'num_sampled', 'source', 'transcript_id'
- output includes all sampled genes in all Chlamydomonas studies from Summary_of_primary_data.csv and omics studies

`05_Incorporate_pathway_information.ipynb`
output: `../../data/intermediate_data_02/sampled_genes.pk`
output columns: 'num_detected', 'num_manipulated', 'num_sampled', 'source', 'transcript_id', 'annotation_version', 'gene_id', 'gene_symbol', 'pathway_id'
- Query pathway information from Phytozome12.1. Only transcripts with identified pathways are returned in the query results.
- As of Jan 06, 2021, cannot query from Phytozomev12.1 anymore. Added `../../data/intermediate_data_02/query_result.csv` to git

`06_update_sampled_genes...ipynb`
output: `../../data/intermediate_data_02/sampled_genes.csv`
output columns: 'num_detected', 'num_manipulated', 'num_sampled', 'source', 'transcript_id', 'annotation_version', 'gene_id', 'gene_symbol', 'pathway_id', 'transcript_id_v5.3.1'
- Add a column of v5.3.1 transcript id to `../../data/intermediate_data_02/sampled_genes.pk`
# Extract relevant genes and generate files in intermediate_data_01/

## Workflow description

1. Run **I. Essentials** and **II. Custom Functions**
2. For each dataset in `lipid_selection/data/raw_data/source_data/`:

    1. Extract basic information:
        - genome_version 
        - database_source 
        - inclusion_criteria 
        - first_author
        - publication_year
      
    2. Append basic information to `lipid_selection/data/intermediate_data_01/basic_info.txt`
        - Use `append_basic_info()`
    
    3. Extract candidate and non-candidate genes.
        - Use `check_excel_data()`, `import_messy_excel()`
    
    4. Export candidate and non-candidate genes to `<first_author>_<year>.txt`.
        - Use `append_genes()` or `append_proteins
        
## Important notes:

* Omit 'Blaby_2013_DS8.xlsx' because it is RNA expression data not transcript expression data.
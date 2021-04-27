## **Files in this folder are supplemental analyses.**

# Summary of files:
    
- **Is the proportion of matches related to year of publication?** 
> Yes, older datasets have more mismatches because the genome version have changed more since then. 

## `01_Intermediate_data_01_exploration.ipynb`
- This notebook only looked at non-selective studies.
**New candidate_gene inclusion criteria:**
- showed up as `True` in at least **one comparison per paper** AND in at least **two papers**

1. **What is the distribution of the number of times each gene is sampled?** 
> Most genes are sampled once or twice.
2. **What is the distribution of the proportion of papers that found >2-fold difference for each gene?** 
> After filtering out genes sampled by only one paper, most genes were detected as candidate genes in one out of eight datasets. 590 genes were detected as candidate genes in at least two of the eight datasets.
3. **Are we looking at substantially diminishing returns of new unique candidate genes per dataset processed already?**
> Randomly sample of all datasets without replacement (n=100 times) shows that the number of unique candidate genes starts to plateau at seven datasets. 
4. **What proportion of TAG synthesis pathway genes are included in this set of unique genes?**
> Four of the 26 TAG synthesis pathway genes are detected as candidate genes.

## `02_Aggregate_experimental_studies_data.ipynb`

Define **selective studies** as studies that 
- Have at these one of the following tags:
> ['knock_out','gene_silence','RT_qPCR','gene_knockout','genetic_transformation']
- have at least 1 gene manipulated.
- are not one of the studies in `match_summary.csv`

*Notice that the difference between `selective_studies.pk` and `selective_studies_reduced.pk` is that `s_source` is a list of `'Firstauthor_year'` in `selective_studies_reduced.pk`. 

Output file: `lipid_selection/data/intermediate_data_02/selective_studies.pk`
Output format:

| transcript_id | s_num_detected | s_num_sampled |s_num_no_effect | s_proportion | s_source                 | annotation_version |
|---------------|------------------|-----------------|------------------|----------------|----------------------------|--------------------|
| `<string>`    | `<int>`          | `<int>`         | `<float>`        | `<float>`      | `'Firstauthor_year'` | `v5.5`             |


Output file: `lipid_selection/data/intermediate_data_02/selective_studies_reduced.pk`
Output format:

| transcript_id | s_num_detected | s_num_sampled |s_num_no_effect | s_proportion | s_source                 | annotation_version |
|---------------|------------------|-----------------|------------------|----------------|----------------------------|--------------------|
| `<string>`    | `<int>`          | `<int>`         | `<float>`        | `<float>`      | `['Firstauthor_year',...]` | `v5.5`             |

**`s_proportion`:** `<float>` 
> `s_proportion` = `s_num_detected`/ `num_selective_source`

**`num_selective_source`:** `<int>` number of selective studies

### Notes:
- Most sampled genes were sampled in one paper only. 
- Most sampled genes have demonstrated significant effect of lipid content in one paper.

## `03_omics_data_collection_simulations.ipynb`
Generates the cumulative unique gene plot for various conditions determining which unique genes to include. 



# DGN: Disease-associated gene network
DGN is a framework integrating scRNA-seq and GWAS data to explore disease-associated genes, cell types, and gene network modules from a network perspective.

## Installation
DGN requires Python >= 3.6 and Java SE Runtime Environment (JRE) >= 1.8. 
```shell
# download 
git clone https://github.com/pmglab/DGN
# install dependent packages
cd DGN
python -m pip install -r requirements.txt
# test
python dgn.py --version
```
## Quick tutorials
### Preparation of DGN Input Files
#### 1. Gene expression matrix of scRNA-seq
For single-cell RNA sequencing (scRNA-seq), the expression matrix must be preprocessed by annotating cell types and then splitting the matrix according to these cell types. 
Each cell type should have a corresponding expression matrix file, where each row represents a gene and each column represents a cell. 
The values should be tab-separated. The file should be named in the format `{CellTypeName}.txt.gz`. All expression matrix files should be placed in a single directory, 
referred to as `EXPR_DIR` in the following demonstrations.

The human and mouse reference single-cell datasets used in our DGN study can be downloaded from https://doi.org/10.6084/m9.figshare.26771971.
#### 2. GWAS Summary Statistics
The GWAS summary statistics should include three columns: chromosome number, position, and p-value. Each row represents a variant. 
The GWAS summary statistics file is referred to as `GWAS_FILE` for demonstration purposes in the following sections.
#### 3. Reference Genotype Data for LD Calculation
Reference genotype data that matches your GWAS population and coordinate version can be downloaded from https://doi.org/10.6084/m9.figshare.26771755. 
For example, the hg19 version of the European population genotype file is named `EUR.hg19.vcf.gz` and will be referred to as `REF_VCF` in the following demonstrations.

### Running DGN
Assume your output directory is set to `OUTPUT_DIR` and the abbreviation for the GWAS phenotype is `PHENO_ABBR`. The command to run is as follows:
```shell
python dgn.py \
--expression_dir EXPR_DIR \
--gwas_path GWAS_FILE \
--vcf_ref REF_VCF \
--output_dir OUTPUT_DIR \
--phenotype_abbr PHENO_ABBR
```
Make sure to replace these placeholder variables with the actual values relevant to your specific situation.
### DGN Output
The output of DGN consists of three main components: disease-associated genes, cell types, and gene network modules.
#### 1. Disease-Associated Genes
#### 2. Disease-Associated Cell Types
#### 3. Disease-Associated Gene Network Modules

## Complete List of DGN Parameters
### Input 1: Gene expression profile
| Flag              | Description  | Default |
|:------------------|:--------------------------|:--------|
| `--expression_dir` | Input a directory that contains gene expression profile of different cell types. Each file represents a specific cell type, with filenames following the convention "{CellTypeName}.txt.gz". |  |
### Input 2: GWAS summary statistics
| Flag            | Description                                                                                   | Default |
|:----------------|:----------------------------------------------------------------------------------------------|:--------|
| `--gwas_path`   | GWAS summary table with header.                                                               |     |
| `--chrom_col`   | Chromosome column in GWAS summary table.                                                      | `CHR`   |
| `--pos_col`     | Base position column in GWAS summary table.                                                   | `BP`    |
| `--p_col`       | P-value column in GWAS summary table.                                                         | `P`     |
| `--buildver`    | Specifies the reference genome version of the coordinates.                                    | `hg19`  |
### Output
| Flag            | Description                                                                                   | Default |
|:----------------|:----------------------------------------------------------------------------------------------|:--------|
| `--output_dir`       | Output directory. |     |
| `--phenotype_abbr`       | Phenotype abbreviation, used as the file identifier for this phenotype in the output results. |      |
### Step 1: Construction of gene co-expression network
| Flag                   | Description                                                                                          | Default   |
|:-----------------------|:-----------------------------------------------------------------------------------------------------|:----------|
| `--trans_gene_symbol`  | Convert the gene ID in the expression profiles to gene symbol of HGNC.                               |       |
| `--min_expr_value`     | The minimum average gene expression level used for constructing co-expression networks.              | `0`       |
| `--min_cell`           | The minimum number of cells or samples required in the analysis.                                     | `100`     |
| `--max_cell`           | The max number of cells or samples required in the analysis.                                         | `1000`    |
| `--min_k`              | The minimum value of k in normalization for co-expression network.                                   | `0.5`     |
| `--max_k`              | The max value of k in normalization for co-expression network.                                       | `1.5`     |
| `--edge_method`        | The method for calculating edge weights (i.e., gene correlations),`pearson` or `cs-core`.            | `cs-core` |
| `--resample_size`              | Sample size for resampling edge weights when correcting the co-expression network.                   | `100000`  |
| `--keep_expr`              | No need to recalculate the gene centrality in the gene co-expression network for the next analysis.  |       |
### Step 2: Disease-associated genes and cell types
| Flag                   | Description                                                                                         | Default   |
|:-----------------------|:----------------------------------------------------------------------------------------------------|:----------|
| `--vcf_ref`  | Specifies a VCF file of genotypes sampled from a reference population.                              |       |
| `--multiple_testing`     | Specifies the method for multiple testing correction. `bonf` denotes performing Bonferroni correction; `benfdr` denotes controlling false discovery rate by the Benjamini–Hochberg method; `fixed` denotes no correction.          | `benfdr`  |
| `--p_value_cutoff`           | Specifies the threshold of the adjusted p-value for fine-mapping.                                    | `0.05`    |
| `--top_n_gene`           | Maximum number of genes entering conditional association analysis.                                        | `1000`    |
| `--nt`              | Specifies the number of threads.                                  | `1`       |
| `--rm_hla`              | Remove HLA region.                                      |       |
| `--edge_method`        | The method for calculating edge weights (i.e., gene correlations),`pearson` or `cs-core`.           | `cs-core` |
| `--java_path`              | Java path.                  | `java`    |
| `--jvm_gb`              | JVM value (GB). | `20`      |

### Step 3: Disease-associated gene network modules
| Flag                   | Description                                                                                                                       | Default |
|:-----------------------|:----------------------------------------------------------------------------------------------------------------------------------|:--------|
| `--assoc_cell_p_cutoff`  | The adjusted p cutoff for associating cell types. The significant cell types will be used for associated network module analysis. | `0.05`  |
| `--module_gene_score_n_top_genes`     | In the module detection function, specify the number of top genes with high disease-related scores for module detection           | `5000`  |
| `--module_cut_edge_weight`           | Minimum edge weight for pruning nodes in modules.                                                                                 | `0.5`   |
| `--module_plot_cut_edge_weight`           | Minimum edge weight for pruning nodes in plotting modules.                                                                        | `0.6`   |
| `--show_node_label_top_n`              | Percentile threshold of degree for displaying the node labels.                                                                    | `5`     |
| `--function_enrich_p_cutoff`              | P-value threshold for functional enrichment analysis by g:Profiler.                                                               | `3`     |
| `--function_enrich_top_n`        | Maximum number of functional enrichment terms displayed.                                        | `3`     |

### Extra analysis: DESE, inferring disease-associated genes and cell types using gene expression.  
| Flag                                | Description                                                                                                                                | Default |
|:------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------|:--------|
| `--run_expr_dese`                   | Calculate the mean expression profile of genes and run DESE.                                                                               | `0.05`  |
| `--normalize_expr`                  | The normalization method for the expression profiles. `no` denote skip normalization, `cpm` denotes count per million (CPM) normalization. | `cpm`   |
| `--log_trans_expr`          | Transform the expression value into log2(x+1).                                                                                         |         |

## Citation





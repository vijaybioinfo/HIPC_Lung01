Data from single-cell T-cell receptor (TCR) repertoire followed by sequencing.
===========

This process was applied to each of several dataset: **(i)** one for lung CD3<sup>+</sup> T cells, **(ii)** one for LLN CD3<sup>+</sup> T cells and **(iii)** one for antigen-reactive CD4<sup>+</sup> T cells from the ARTE assay. This section comprises the following major steps:

* Preprocessing of reads by mapping them to the reference genome and TCR clonotype assembly in a library-wise manner (with `cellranger`, v7.1).
* Aggregation of TCR clonotypes across libraries for selected groups: **(i)** one for lung CD3<sup>+</sup> T cells, **(ii)** one for LLN CD3<sup>+</sup> T cells and **(iii)** one for antigen-reactive CD4<sup>+</sup> T cells from the ARTE assay.
* Clonotype-wise summarization: encompassing CDR3 and V and J gene segments, absolute and relative frequencies in a subset-wise fashion and across the board, etc.

The first step is described in the methods section of our manuscript in detail. Here we provide the code for the rest of the steps.

---
# Scripts and job scripts
The following scripts are located under the `./gen_scripts` folder, extensively self-documented and with a brief description below. The aggregation-specific job scripts ran (and their corresponding standard outputs) can be found in the `aggr_vdj` folder under each aggregation folder.
* `./gen_scripts/aggr_vdj.2.2.R`: Performs aggregation of TCR clonotypes across libraries for selected groups, as defined above.
* `./gen_scripts/comp_tcr_info.0.1.R`: To summarize clonotype-wise information, also as described above. This was applied specifically to the lung and LLN CD3<sup>+</sup> T cell aggregations only. This was not applied to the aggregation of antigen-reactive CD4<sup>+</sup> TCRs; instead, we proceeded with a dataset-specific script provided in the main `./../../reference_set_def` folder.

---
# System requirements and installation
* R (version $\geq$ v3.6).
* `Seurat` v3.2.2
* `data.table`
* `stringr`
* `tidyr`

---
# Data requirements
For the second step (aggregation), the outputs from `cellranger vdj` (v7.1) are required, as described in detail elsewhere. To obtain the raw data for the results presented in our study, you must download it from our GEO submission (to be provided).<br>
For the third step (summarization), the outputs from the aggregation, as well as a Seurat object are needed. The latter is the output from the preprocessing of the scRNA-seq (transcriptome) data, as described in `./../transcriptome_data`.
# Evolution of T cell responses as a function of multiple COVID-19 boosters
------------

Description
------------

The long-term effects of repeated COVID-19 vaccinations on adaptive immunity remain incompletely understood. Here, we conducted a comprehensive three-year longitudinal study examining T cell and antibody responses in 68 vaccinated individuals without reported symptomatic infections. We observed distinct dynamics in humoral and cellular immune responses across multiple vaccine doses. While antibody titers incrementally increased and stabilized with each booster, T cell responses rapidly plateaued, maintaining remarkable stability across CD4+ and CD8+ subsets. Notably, approximately 30% of participants showed CD4+ T cell reactivity to non-Spike antigens, consistent with asymptomatic infections. Single-cell RNA sequencing revealed a diverse landscape of Spike-specific T cell phenotypes, with no evidence of increased exhaustion or significant functional impairment. However, qualitative changes were observed in individuals with evidence of asymptomatic infection, exhibiting unique immunological characteristics, including increased frequencies of Th17-like CD4+ T cells and GZMKhi/IFNR-like CD8+ T cell subsets. Remarkably, repeated vaccinations in this group were associated with a progressive increase in regulatory T cells, potentially indicating a balanced immune response that may mitigate immunopathology. By stimulating T cell memory, boosters contribute to a stable and enhanced immune response, which may provide better protection against symptomatic infections.

This repository contains the data and scripts used to analyze the scRNA-seq data produced in our study and produce some of the associated figures shown on our manuscript.

Requirements
------------

This project was done using the following modules/programs:

* [R](https://cran.r-project.org/) (v3.6.1)
* [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) (v7.1.0)
* [Seurat](https://satijalab.org/seurat) (v3.2.3)
* [MAST](https://github.com/RGLab/MAST) (v1.12.0)

Raw data
------------
* The single-cell RNA-seq raw and processed files from our study can be downloaded via GEO with accession number to be provided soon.

Raw data processing  
------------

### Single-cell
* Preprocesing was performed with 10x's cellranger (v7.1.0) as generally described (see Methods section of publication).
* Quality control (QC), unbiased clustering, dimensionality reduction and cluster annotation were performed with the workflow detailed in [this repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3)

> Relevant scripts are located in: ./pre-processing  

For more specific information about the data generation and processing, please check the Methods section of our manuscript.  

Figures
------------
> Relevant scripts are located in: ./figures


Downstream Analysis
------------
* DGEA - You can follow [MAST's docs](https://github.com/RGLab/MAST)
> Relevant scripts all located in: ./downstream_analysis


Usage & Citation
--------------

If you want to clone this repository run:
```bash
git clone https://github.com/vijaybioinfo/PBT_2023.git
```
Please cite the following manuscript if you are using this repository:  
To update upon publication.

Maintainers
-----------

Current maintainer:
* Vicente Fajardo Rosas (vfajardo@lji.org) 

Sette & Vijayanand Labs.  
Center for Infectious Disease and Vaccine Research, La Jolla Institute for Immunology La Jolla, CA 92037, USA


Contact
-----------
Please email any of the current maintainers.

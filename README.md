# The human lung is a reservoir of tissue-resident memory T cells against a broad spectrum of pathogens
------------

Description
------------

Repo for our manuscript titled as above. This repository self-contains some and links to other repos with further scripts used to process and analyze the data produced in our study. We also provide the code we wrote to produce the associated figures and tables.

### Summary
Immune responses in the lung are essential for preventing and controlling a wide range of infections. Tissue-resident memory T (TRM) cells are critical for frontline immunity—yet in the murine lung they undergo rapid attrition to preserve gas exchange, accumulating in the draining lymph nodes (LLN) for resupply through retrograde migration. Whether this paradigm applies to humans is currently unknown. Here, we present the first comprehensive analysis of T cells from human lung and its LLN, defining the prevalence and properties of lung TRM cells specific to a broad spectrum of pathogens. Using a T-cell receptor (TCR)-guided approach that integrates single-cell transcriptomics with paired TCR repertoire profiling, we mapped the pathogen-specificity of over 87,000 lung T cells from 40 individuals, the majority of whom harbored TRM cell populations specific to multiple viruses. Thus, in contrast to mice, the human lung itself serves as a reservoir of clonally-expanded TRM cells specific to a broad spectrum of pathogens. Although the LLN contained a small population of TRM cells, over 80% of the TCRs from highly expanded TRM clones in the lung were absent in the LLN, suggesting that they do not serve as the major reservoir for lung TRM cells. We also found that the prevalence and properties of pathogen-specific lung TRM cells are influenced more strongly by donor-intrinsic factors than pathogen identity. This suggests that the design of effective vaccines targeting airborne pathogens needs to account for host-intrinsic variability in the generation and maintenance of tissue-resident memory cells. Overall, our findings reveal that the human lung retains a stable and diverse pool of pathogen-specific TRM cells, suggesting that strategies to bolster these responses could provide durable protection against severe lung infections.


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
git clone https://github.com/vijaybioinfo/HIPC_Lung01.git
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

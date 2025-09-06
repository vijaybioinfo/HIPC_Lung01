Data from hashtag oligonucleotide (HTO) by feature barcoding, followed by sequencing.
===========

This process was applied to each HTO sample, to assign donor identities to a large fraction of all cell barcodes and identify potential multiplets. This comprises the following major steps:

* Preprocessing of reads by mapping them to the reference feature barcoding reference in a library (sample)-wise manner (with `cellranger`, v7.1).
* Assignment of donor IDs to cell barcodes ("donor deconvolution").

The first step is described in the methods section of our manuscript in detail. Here we provide the code for the other step.

---
# Scripts and job scripts
The following scripts are located under the `./gen_scripts` folder, extensively self-documented and with a brief description below. The library-specific job scripts ran (and their corresponding standard outputs) can be found in the `./jobs_scripts` folder under this repo folder.
* `./gen_scripts/HTO_based_demultiplexing.2.4.R`: Performs donor deconvolution, as defined above and explained in more detail elsewhere<sup>[1](https://link.springer.com/article/10.1186/s13059-018-1603-1),[2](https://www.nature.com/articles/s41592-019-0433-8)</sup>

---
# System requirements and installation
* R (version $\geq$ v3.6).
* `Seurat` v3.2.2
* `stringr`
* `ggplot2`

---
# Data requirements
The outputs from `cellranger count` adapted for feature barcoding (v7.1) are required. To obtain the raw data for the results presented in our study, you must download it from our GEO submission (to be provided).<br>
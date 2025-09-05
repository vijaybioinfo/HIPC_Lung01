Data preprocessing
===========

Description
------------

Sequencing reads from scRNA-seq, scTCR-seq and feature barcoding (hashtag oligonucleotide (HTO)) data were demultiplexed and aligned to references using Cell Ranger (v7.1.0, 10x Genomics). See the Methods sections of our manuscript for full details.<br/>

Details on the preprocessing of each different data modality are provided in the separate folders, namely for:
* Single-cell RNA-seq data (`transcriptome_data`).
* Single-cell feature barcoding (hashtag oligonucleotide, HTO) data (`hto_data`).
* Single-cell T-cell receptor (TCR) repertoire data (`sc_tcr_data`).
* Bulk TCR data (`bulk_tcr_data`).
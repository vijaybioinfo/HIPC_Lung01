#!/bin/sh
#SBATCH --job-name=BOOST_Round-4_CD8_SA
#SBATCH --output=/mnt/hpcscratch/vfajardo/BOOST/jobs_scripts/seurat_analysis//BOOST_CD8/BOOST_Round-4_CD8/BOOST_Round-4_CD8_09-06-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag/BOOST_Round-4_CD8_general_seurat_analysis.out.txt
#SBATCH --error=/mnt/hpcscratch/vfajardo/BOOST/jobs_scripts/seurat_analysis//BOOST_CD8/BOOST_Round-4_CD8/BOOST_Round-4_CD8_09-06-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag/BOOST_Round-4_CD8_general_seurat_analysis.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem=380g
#SBATCH --time=40:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: BOOST_Round-4_CD8"
echo -e "Job output directory: /mnt/hpcscratch/vfajardo/BOOST/jobs_scripts/seurat_analysis//BOOST_CD8/BOOST_Round-4_CD8/BOOST_Round-4_CD8_09-06-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Seurat script version (absolute path): /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.3.9.R"
echo -e "Project name: BOOST_Round-4_CD8"
echo -e "Absolute path to count matrix file: "
echo -e "Feature ID: name"
echo -e "Minimum and maximum thresholds for counts number per cell: 1500 and 20000"
echo -e "Minimum and maximum thresholds for genes number per cell: 800 and 4400"
echo -e "Maximum mitochondrial genes percentage per cell: 10"
echo -e "Normalization method: LogNormalize"
echo -e "Variable features selection method: vst"
echo -e "Number of most variable genes to be considered: 30"
echo -e "Resolution argument for cluster analysis: 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)'"
echo -e "PCs to be taken into account for dimensionality reduction: 'c(30)'\n\n"
echo -e "Mean cutoff: 0.01"
echo -e "Should likely low quality cells be filtered out? TRUE"
echo -e "Absolute path to markers file (or NULL value): /home/vfajardo/scripts/seurat_analysis/general_data/CD4_Tcell_markers.1.3.1.RData"
echo -e "### --------------------- Annotations arguments --------------------- ###"
echo -e "Should annotations be added? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Absolute path to annotations table: /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/08-07-2023/aggr/data/BOOST_Round-4_CD8_aggr_table_annotated.0.2.csv"
echo -e "Lane ID: 'cell.type.tag;facs.sorting.batch.tag'"
echo -e "10X Chromium Batch ID: chrom.batch.tag"
echo -e "Sequencing batch ID: seq.batch.tag"
echo -e "HTO ID: hashtag.tag"
echo -e "Demuxlet ID: "
echo -e "Overlay ID: donor.tag"
echo -e "### ---------------------- Filtering arguments ---------------------- ###"
echo -e "Should there be any kind of preprocess filetring? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Tags-related criteria: --TagsCriteria /mnt/bioadhoc-temp/Groups/vd-sette/BOOST/paper_developments/BOOST_1/data_for_subsetting/subset_hto-all_pp-S_cohort-all/TagsSubsetCriteria.csv"
echo -e "Features-related criteria: "

echo -e "Here starts analysis with R package seurat...\n\n"

module load R/3.6.1

Rscript /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.3.9.R --ReportsPath /mnt/hpcscratch/vfajardo/BOOST/seurat_analysis//BOOST_CD8/BOOST_Round-4_CD8/BOOST_Round-4_CD8_09-06-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag --PrjName BOOST_Round-4_CD8 --DataFile /mnt/hpcscratch/vfajardo/sequencing_data/08-07-2023/aggr/BOOST_Round-4_CD8/outs/count/filtered_feature_bc_matrix --Raw10x /mnt/hpcscratch/vfajardo/sequencing_data/08-07-2023/aggr/BOOST_Round-4_CD8/outs/count/filtered_feature_bc_matrix --InputType matrix   --FeatureID name --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 10 --GenNormMethod LogTransform --NormMethod LogNormalize --FVFsMethod vst --FeatsForDSA 30 --PCs 'c(30)' --ForHarmony 'c("seq.batch.tag")' --Resolution 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)' --IsFivePrime FALSE --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")'  --FilterOut TRUE --MarkersFile /home/vfajardo/scripts/seurat_analysis/general_data/CD4_Tcell_markers.1.3.1.RData --DEA FALSE --DimReduction 'c("umap")' --RAM 380 --DoAnnotate TRUE --AggrTable /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/08-07-2023/aggr/data/BOOST_Round-4_CD8_aggr_table_annotated.0.2.csv --DonorsMetaData /mnt/bioadhoc-temp/Groups/vd-sette/BOOST/paper_developments/BOOST_1/donors_metadata/DonorMetadata_ProcessedTable_Round-4.csv  --MergeBy donor.tag --LaneID 'cell.type.tag;facs.sorting.batch.tag' --MultAnnPerLane TRUE --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag  --DoSubset TRUE --TagsCriteria /mnt/bioadhoc-temp/Groups/vd-sette/BOOST/paper_developments/BOOST_1/data_for_subsetting/subset_hto-all_pp-S_cohort-all/TagsSubsetCriteria.csv 

echo -e "Job completed!\nCheck for errors if any."

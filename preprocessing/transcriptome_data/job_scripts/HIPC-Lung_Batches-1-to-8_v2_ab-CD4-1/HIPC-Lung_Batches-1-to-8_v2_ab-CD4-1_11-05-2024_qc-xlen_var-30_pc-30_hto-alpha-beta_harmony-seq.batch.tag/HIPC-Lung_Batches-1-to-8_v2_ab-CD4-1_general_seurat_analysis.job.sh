#!/bin/sh
#SBATCH --job-name=HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_SA
#SBATCH --output=/path/to/user/HIPC/jobs_scripts/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_general_seurat_analysis.out.txt
#SBATCH --error=/path/to/user/HIPC/jobs_scripts/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_general_seurat_analysis.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=400g
#SBATCH --time=120:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1"
echo -e "Job output directory: /path/to/user/HIPC/jobs_scripts/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Seurat script version (absolute path): /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.4.3.R"
echo -e "Project name: HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1"
echo -e "Absolute path to count matrix file: "
echo -e "Feature ID: name"
echo -e "Minimum and maximum thresholds for counts number per cell: 500 and 10000"
echo -e "Minimum and maximum thresholds for genes number per cell: 400 and 3000"
echo -e "Maximum mitochondrial genes percentage per cell: 15"
echo -e "Normalization method: LogNormalize"
echo -e "Variable features selection method: vst"
echo -e "Number of most variable genes to be considered: 30"
echo -e "Resolution argument for cluster analysis: 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)'"
echo -e "PCs to be taken into account for dimensionality reduction: 'c(30)'\n\n"
echo -e "Mean cutoff: 0.01"
echo -e "Should likely low quality cells be filtered out? TRUE"
echo -e "Absolute path to markers file (or NULL value): /home/vfajardo/scripts/seurat_analysis/general_data/T-cell_markers.1.0.RData"
echo -e "### --------------------- Annotations arguments --------------------- ###"
echo -e "Should annotations be added? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Absolute path to annotations table: /path/to/user/sequencing_data/11-17-2023/aggr/data/HIPC-Lung_Batches-1-to-8_aggr_table_annotated.0.2.csv"
echo -e "Lane ID: 'library.id.tag;pre.donor.id.tag;facs.sorting.batch.tag'"
echo -e "10X Chromium Batch ID: chrom.batch.tag"
echo -e "Sequencing batch ID: seq.batch.tag"
echo -e "HTO ID: hashtag.tag"
echo -e "Demuxlet ID: "
echo -e "Overlay ID: donor.tag"
echo -e "### ---------------------- Filtering arguments ---------------------- ###"
echo -e "Should there be any kind of preprocess filetring? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Tags-related criteria: --TagsCriteria /path/to/user/HIPC/paper_developments/HIPC-Lung/data_for_subsetting/subset_hto-alpha-beta/TagsSubsetCriteria.csv"
echo -e "Features-related criteria: "

echo -e "Here starts analysis with R package seurat...\n\n"

module load R/3.6.1

Rscript /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.4.3.R \
--ReportsPath /root/hpcscratch/vfajardo/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag --PrjName HIPC-Lung_Batches-1-to-8_v2_ab-CD4-1 --RAM 400 \
--DataFile /path/to/user/HIPC/seurat_analysis/HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_CD4-3/HIPC-Lung_Batches-1-to-8_v2_CD4-3_11-04-2024_qc-xlen_var-30_pc-30_hto-custom-1_harmony-seq.batch.tag/seurat_objects/SeuratObjectForPrjHIPC-Lung_Batches-1-to-8_v2_CD4-3_WithArgs_NoPCs_30.RDS  --InputType seurat --FeatureID name --DoPreSubset TRUE --PreSubsetCriteria /path/to/user/HIPC/paper_developments/HIPC-Lung/data_for_subsetting/custom_subsets//HIPC-Lung//HIPC-Lung_Batches-1-to-8_v2_CD4-3/HIPC-Lung_Batches-1-to-8_v2_CD4-3_11-04-2024_qc-xlen_var-30_pc-30_hto-custom-1_harmony-seq.batch.tag/ab-CD4-1/TagsSubsetCriteria.csv \
--FilterOut TRUE --minCounts 500 --maxCounts 10000 --minFeatures 400 --maxFeatures 3000 --maxMP 15 --IsFivePrime FALSE \
--GenNormMethod LogTransform --NormMethod LogNormalize \
--FVFsMethod vst --FeatsForDSA 30 --MeanCutoff 0.01 --PropCutoff 0.001 \
--PCs 'c(30)' --ForHarmony 'c("seq.batch.tag")' --VarsToRegress 'c("nCount_RNA", "percent.mt")'  --DimReduction 'c("umap")' --Resolution 'c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)' \
--MarkersFile /home/vfajardo/scripts/seurat_analysis/general_data/T-cell_markers.1.0.RData --DEA FALSE \
--DoAnnotate TRUE --AggrTable /path/to/user/sequencing_data/11-17-2023/aggr/data/HIPC-Lung_Batches-1-to-8_aggr_table_annotated.0.2.csv   --MergeBy donor.tag --LaneID 'library.id.tag;pre.donor.id.tag;facs.sorting.batch.tag' --MultAnnPerLane TRUE --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag  --NewMeta /path/to/user/HIPC/paper_developments/HIPC-Lung/exploratory_analyses/state_deconvolution/state_deconvolution_2024-10-30/final_assignments/HIPC-40_ExVivo_FinalAssignments.csv \
--DoSubset TRUE --TagsCriteria /path/to/user/HIPC/paper_developments/HIPC-Lung/data_for_subsetting/subset_hto-alpha-beta/TagsSubsetCriteria.csv 

echo -e "Job completed!\nCheck for errors if any."

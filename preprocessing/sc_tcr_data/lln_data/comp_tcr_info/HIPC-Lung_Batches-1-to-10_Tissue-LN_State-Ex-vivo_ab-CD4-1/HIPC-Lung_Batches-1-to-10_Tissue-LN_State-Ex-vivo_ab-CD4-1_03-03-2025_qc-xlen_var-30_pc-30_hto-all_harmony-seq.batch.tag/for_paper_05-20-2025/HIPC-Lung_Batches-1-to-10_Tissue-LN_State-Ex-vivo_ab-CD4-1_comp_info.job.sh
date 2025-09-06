#!/bin/sh
#SBATCH --job-name=HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_TCR
#SBATCH --output=/path/to/user/HIPC/jobs_scripts/TCR_data_analysis/comp_tcr_info//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/for_paper_05-20-2025/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_comp_info.out.txt
#SBATCH --error=/path/to/user/HIPC/jobs_scripts/TCR_data_analysis/comp_tcr_info//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/for_paper_05-20-2025/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_comp_info.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=150g
#SBATCH --time=120:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"

echo -e "Here starts analysis with R...\n\n"

module load R/3.6.1

Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/get_comp_tcr_info.0.1.R --ReportsPath /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/comp_tcr_info/for_paper --AggrPath /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025 --GExData /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/seurat_objects/SeuratObjectForPrjHIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_WithArgs_NoPCs_30.RDS  --ClustsLab custom.cluster.tag --VDJAggrTable /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025/vdj_aggr_table.csv --GExAggrTable /path/to/user/sequencing_data/12-15-2024/aggr/data/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_aggr_table_annotated.0.1.csv  --FeatureID name --FiltCriteria /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/comp_tcr_info/for_paper/input_data/TagsCriteria.csv --TransformRules /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_ab-CD4-1_03-03-2025_qc-xlen_var-30_pc-30_hto-all_harmony-seq.batch.tag/comp_tcr_info/for_paper/input_data/RulesForTagTransformation.csv


echo -e "Job completed!\nCheck for errors if any."

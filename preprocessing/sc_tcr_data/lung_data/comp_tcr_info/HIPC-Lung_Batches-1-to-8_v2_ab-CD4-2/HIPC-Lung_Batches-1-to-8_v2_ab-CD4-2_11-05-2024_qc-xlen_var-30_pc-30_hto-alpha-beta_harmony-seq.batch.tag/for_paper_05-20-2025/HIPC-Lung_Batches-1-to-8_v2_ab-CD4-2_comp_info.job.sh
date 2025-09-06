#!/bin/sh
#SBATCH --job-name=HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_TCR
#SBATCH --output=/path/to/user/HIPC/jobs_scripts/TCR_data_analysis/comp_tcr_info//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/for_paper_05-20-2025/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_comp_info.out.txt
#SBATCH --error=/path/to/user/HIPC/jobs_scripts/TCR_data_analysis/comp_tcr_info//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/for_paper_05-20-2025/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_comp_info.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=150g
#SBATCH --time=120:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"

echo -e "Here starts analysis with R...\n\n"

module load R/3.6.1

Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/get_comp_tcr_info.0.1.R --ReportsPath /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/comp_tcr_info/for_paper --AggrPath /path/to/user/sequencing_data/11-17-2023/aggr_vdj/HIPC-Lung_Batches-1-to-8_v2/HIPC-Lung_Batches-1-to-8_v2_1.0_02-12-2024 --GExData /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/seurat_objects/SeuratObjectForPrjHIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_WithArgs_NoPCs_30.RDS --ClustsLab RNA_snn_res.0.2 --VDJAggrTable /path/to/user/sequencing_data/11-17-2023/aggr_vdj/HIPC-Lung_Batches-1-to-8_v2/HIPC-Lung_Batches-1-to-8_v2_1.0_02-12-2024/vdj_aggr_table.csv --GExAggrTable /path/to/user/sequencing_data/11-17-2023/aggr/data/HIPC-Lung_Batches-1-to-8_v2_aggr_table_annotated.0.2.csv --FiltCriteria /path/to/user/HIPC/seurat_analysis//HIPC-Lung/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2/HIPC-Lung_Batches-1-to-8_v2_ab-CD4-2_11-05-2024_qc-xlen_var-30_pc-30_hto-alpha-beta_harmony-seq.batch.tag/comp_tcr_info/for_paper/input_data/TagsCriteria.csv


echo -e "Job completed!\nCheck for errors if any."

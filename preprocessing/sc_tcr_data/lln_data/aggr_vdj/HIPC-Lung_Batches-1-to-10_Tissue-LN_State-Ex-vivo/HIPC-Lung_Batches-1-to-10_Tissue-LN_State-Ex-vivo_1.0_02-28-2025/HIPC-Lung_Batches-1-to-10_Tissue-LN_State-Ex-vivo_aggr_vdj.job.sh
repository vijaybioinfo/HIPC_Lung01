#!/bin/sh
#SBATCH --job-name=_aggr_vdj
#SBATCH --output=/path/to/user/sequencing_data/12-15-2024/jobs_scripts/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_aggr_vdj.out.txt
#SBATCH --error=/path/to/user/sequencing_data/12-15-2024/jobs_scripts/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_aggr_vdj.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200g
#SBATCH --time=40:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo"
echo -e "Job output directory: /path/to/user/sequencing_data/12-15-2024/jobs_scripts/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Reports: /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025"
echo -e "General input-data path: /path/to/user/sequencing_data/12-15-2024/vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo"
echo -e "Path to table describing samples to be aggregated: /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/data/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_vdj_aggr_table.1.0.csv"

module load R/3.6.1

Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/aggr_vdj.2.2.R --ReportsPath /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_1.0_02-28-2025 --GenInputPath /path/to/user/sequencing_data/12-15-2024/vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo --AggrTable /path/to/user/sequencing_data/12-15-2024/aggr_vdj/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo/data/HIPC-Lung_Batches-1-to-10_Tissue-LN_State-Ex-vivo_vdj_aggr_table.1.0.csv --SampleCountThold 2 --FreqThold 2

echo -e "Job completed!\nCheck for errors if any."

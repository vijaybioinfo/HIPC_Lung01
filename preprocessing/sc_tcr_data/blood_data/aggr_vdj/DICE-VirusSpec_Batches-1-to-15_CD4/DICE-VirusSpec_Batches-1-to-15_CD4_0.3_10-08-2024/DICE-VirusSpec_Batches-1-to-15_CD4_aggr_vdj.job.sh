#!/bin/sh
#SBATCH --job-name=_aggr_vdj
#SBATCH --output=/path/to/user/sequencing_data/08-21-2023/jobs_scripts/aggr_vdj/DICE-VirusSpec_Batches-1-to-15_CD4/DICE-VirusSpec_Batches-1-to-15_CD4_0.3_10-08-2024/DICE-VirusSpec_Batches-1-to-15_CD4_aggr_vdj.out.txt
#SBATCH --error=/path/to/user/sequencing_data/08-21-2023/jobs_scripts/aggr_vdj/DICE-VirusSpec_Batches-1-to-15_CD4/DICE-VirusSpec_Batches-1-to-15_CD4_0.3_10-08-2024/DICE-VirusSpec_Batches-1-to-15_CD4_aggr_vdj.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=180g
#SBATCH --time=40:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: DICE-VirusSpec_Batches-1-to-15_CD4"
echo -e "Job output directory: /path/to/user/sequencing_data/08-21-2023/jobs_scripts/aggr_vdj/DICE-VirusSpec_Batches-1-to-15_CD4/DICE-VirusSpec_Batches-1-to-15_CD4_0.3_10-08-2024"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Reports: /path/to/user/sequencing_data/08-21-2023/aggr_vdj/DICE-VirusSpec_Batches-1-to-15_CD4/DICE-VirusSpec_Batches-1-to-15_CD4_0.3_10-08-2024"
echo -e "General input-data path: /path/to/user/sequencing_data/08-21-2023/vdj/DICE-VirusSpec_Batches-1-to-15_CD4"
echo -e "Path to table describing samples to be aggregated: /path/to/user/sequencing_data/08-21-2023/aggr_vdj/data/DICE-VirusSpec_Batches-1-to-15_CD4_aggr_table.0.3.csv"

module load R/3.6.1

Rscript /home/vfajardo/scripts/TCR_data_analysis/10X_preliminary_analyses/aggr_vdj.2.2.R --ReportsPath /path/to/user/sequencing_data/08-21-2023/aggr_vdj/DICE-VirusSpec_Batches-1-to-15_CD4/DICE-VirusSpec_Batches-1-to-15_CD4_0.3_10-08-2024 --GenInputPath /path/to/user/sequencing_data/08-21-2023/vdj/DICE-VirusSpec_Batches-1-to-15_CD4 --AggrTable /path/to/user/sequencing_data/08-21-2023/aggr_vdj/data/DICE-VirusSpec_Batches-1-to-15_CD4_aggr_table.0.3.csv --SampleCountThold 2 --FreqThold 2

echo -e "Job completed!\nCheck for errors if any."

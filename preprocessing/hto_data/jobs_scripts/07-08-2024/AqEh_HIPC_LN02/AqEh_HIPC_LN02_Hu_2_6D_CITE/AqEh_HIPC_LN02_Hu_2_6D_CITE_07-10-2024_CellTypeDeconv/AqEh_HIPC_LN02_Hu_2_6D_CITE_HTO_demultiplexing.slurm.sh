#!/bin/sh
#SBATCH --job-name=AqEh_HIPC_LN02_Hu_2_6D_CITE_HTODem
#SBATCH --output=/path/to/user/HIPC/jobs_scripts/deconvolution/HTO_based/07-08-2024/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/AqEh_HIPC_LN02_Hu_2_6D_CITE_07-10-2024_CellTypeDeconv/AqEh_HIPC_LN02_Hu_2_6D_CITE_HTO_demultiplexing.out.txt
#SBATCH --error=/path/to/user/HIPC/jobs_scripts/deconvolution/HTO_based/07-08-2024/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/AqEh_HIPC_LN02_Hu_2_6D_CITE_07-10-2024_CellTypeDeconv/AqEh_HIPC_LN02_Hu_2_6D_CITE_HTO_demultiplexing.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30g
#SBATCH --time=02:00:00

echo -e "\n\n######### --------- Job to Run HTO-based deconvolution --------- #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: AqEh_HIPC_LN02_Hu_2_6D_CITE"
echo -e "Job output directory: /path/to/user/HIPC/jobs_scripts/deconvolution/HTO_based/07-08-2024/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/AqEh_HIPC_LN02_Hu_2_6D_CITE_07-10-2024_CellTypeDeconv"
echo -e "### --------------- HTO-based deconvolution arguments --------------- ###"
echo -e "Absolute path to reports directory: /path/to/user/HIPC/deconvolution/HTO_based/07-08-2024/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/AqEh_HIPC_LN02_Hu_2_6D_CITE_07-10-2024_CellTypeDeconv"
echo -e "Absolute path to HTO data: /path/to/user/sequencing_data/07-08-2024/count/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/outs/raw_feature_bc_matrix"
echo -e "Gene expression data or NULL value: "
echo -e "UMI threshold for raw data filtering or NULL value: 100"

module load R/3.6.1

Rscript /home/vfajardo/scripts/deconvolution/HTO_based/HTO_based_demultiplexing.2.4.R --PrjName AqEh_HIPC_LN02_Hu_2_6D_CITE --ReportsPath /path/to/user/HIPC/deconvolution/HTO_based/07-08-2024/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/AqEh_HIPC_LN02_Hu_2_6D_CITE_07-10-2024_CellTypeDeconv --HTOData /path/to/user/sequencing_data/07-08-2024/count/AqEh_HIPC_LN02/AqEh_HIPC_LN02_Hu_2_6D_CITE/outs/raw_feature_bc_matrix --HTONames "c('CD4', 'CD8A')" --FCThold 3  --UMIThold 100 --DoReclass doublet

echo -e "Job completed!\nCheck for errors if any."

#!/bin/sh
#SBATCH --job-name=AqEh06_Hu_TCL11_2D_CITE_HTODem
#SBATCH --output=/root/bioadhoc-temp/Groups/vd-vijay/vfajardo/HIPC/jobs_scripts/deconvolution/HTO_based/06-28-2023/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/AqEh06_Hu_TCL11_2D_CITE_07-25-2023_CellTypeDeconv/AqEh06_Hu_TCL11_2D_CITE_HTO_demultiplexing.out.txt
#SBATCH --error=/root/bioadhoc-temp/Groups/vd-vijay/vfajardo/HIPC/jobs_scripts/deconvolution/HTO_based/06-28-2023/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/AqEh06_Hu_TCL11_2D_CITE_07-25-2023_CellTypeDeconv/AqEh06_Hu_TCL11_2D_CITE_HTO_demultiplexing.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30g
#SBATCH --time=02:00:00

echo -e "\n\n######### --------- Job to Run HTO-based deconvolution --------- #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: AqEh06_Hu_TCL11_2D_CITE"
echo -e "Job output directory: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/HIPC/jobs_scripts/deconvolution/HTO_based/06-28-2023/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/AqEh06_Hu_TCL11_2D_CITE_07-25-2023_CellTypeDeconv"
echo -e "### --------------- HTO-based deconvolution arguments --------------- ###"
echo -e "Absolute path to reports directory: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/HIPC/deconvolution/HTO_based/06-28-2023/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/AqEh06_Hu_TCL11_2D_CITE_07-25-2023_CellTypeDeconv"
echo -e "Absolute path to HTO data: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/06-28-2023/count/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/outs/raw_feature_bc_matrix"
echo -e "Gene expression data or NULL value: "
echo -e "UMI threshold for raw data filtering or NULL value: 100"

module load R/3.6.1

Rscript /home/vfajardo/scripts/deconvolution/HTO_based/HTO_based_demultiplexing.2.4.R --PrjName AqEh06_Hu_TCL11_2D_CITE --ReportsPath /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/HIPC/deconvolution/HTO_based/06-28-2023/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/AqEh06_Hu_TCL11_2D_CITE_07-25-2023_CellTypeDeconv --HTOData /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/06-28-2023/count/AqEh_HIPC06/AqEh06_Hu_TCL11_2D_CITE/outs/raw_feature_bc_matrix --HTONames "c('CD4', 'CD8A')" --FCThold 3  --UMIThold 100 --DoReclass doublet

echo -e "Job completed!\nCheck for errors if any."

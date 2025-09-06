#!/bin/sh
#SBATCH --job-name=DICE_TCR014_Hu_1_14_CITE_HTODem
#SBATCH --output=/root/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/jobs_scripts/deconvolution/HTO_based/04-23-2023/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/DICE_TCR014_Hu_1_14_CITE_05-02-2023_100umis_thold/DICE_TCR014_Hu_1_14_CITE_HTO_demultiplexing.out.txt
#SBATCH --error=/root/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/jobs_scripts/deconvolution/HTO_based/04-23-2023/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/DICE_TCR014_Hu_1_14_CITE_05-02-2023_100umis_thold/DICE_TCR014_Hu_1_14_CITE_HTO_demultiplexing.err.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus=3
#SBATCH --mem=30g
#SBATCH --time=02:00:00

echo -e "\n\n######### --------- Job to Run HTO-based deconvolution --------- #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: DICE_TCR014_Hu_1_14_CITE"
echo -e "Job output directory: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/jobs_scripts/deconvolution/HTO_based/04-23-2023/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/DICE_TCR014_Hu_1_14_CITE_05-02-2023_100umis_thold"
echo -e "### --------------- HTO-based deconvolution arguments --------------- ###"
echo -e "Absolute path to reports directory: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/deconvolution/HTO_based/04-23-2023/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/DICE_TCR014_Hu_1_14_CITE_05-02-2023_100umis_thold"
echo -e "Absolute path to HTO data: /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/04-23-2023/count/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/outs/raw_feature_bc_matrix"
echo -e "Gene expression data or NULL value: "
echo -e "UMI threshold for raw data filtering or NULL value: 100"

Rscript /home/vfajardo/scripts/deconvolution/HTO_based/HTO_based_demultiplexing.2.4.R --PrjName DICE_TCR014_Hu_1_14_CITE --ReportsPath /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/deconvolution/HTO_based/04-23-2023/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/DICE_TCR014_Hu_1_14_CITE_05-02-2023_100umis_thold --HTOData /root/bioadhoc-temp/Groups/vd-vijay/vfajardo/sequencing_data/04-23-2023/count/DICE_TCR014/DICE_TCR014_Hu_1_14_CITE/outs/raw_feature_bc_matrix --HTONames NULL --FCThold 3  --UMIThold 100 --DoReclass doublet

echo -e "Job completed!\nCheck for errors if any."

#!/bin/bash
#SBATCH --job-name=DTB15_v2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200g
#SBATCH --time=168:00:00
#SBATCH --output=/path/to/user/R24/paper_developments/DICE-VirusSpec/exploratory_analyses/preliminary_process/jobs_scripts/run_prelim_process_v0.3_2024-10-14/run_prelim_process.out.txt
#SBATCH --error=/path/to/user/R24/paper_developments/DICE-VirusSpec/exploratory_analyses/preliminary_process/jobs_scripts/run_prelim_process_v0.3_2024-10-14/run_prelim_process.err.txt

module unload R
module load R/3.6.1

Rscript /path/to/user/R24/paper_developments/DICE-VirusSpec/exploratory_analyses/preliminary_process/jobs_scripts/prelim_process.3.0.R --BatchesLab 'Batches-1-to-15' --SeqDate '08-21-2023' --AggrVers '0.3' --AggrDate '10-08-2024'
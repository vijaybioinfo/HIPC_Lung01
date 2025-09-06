#!/bin/bash
#SBATCH --job-name=PAPP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=150g
#SBATCH --time=160:00:00
#SBATCH --output=/path/to/user/HIPC/paper_developments/HIPC-Lung/exploratory_analyses/prol_assay_reps/jobs_scripts/jobs/run_prol_assay_reps_art-thold-0.5_2024-12-31/run_prol_assay_tcr_process.out.txt
#SBATCH --error=/path/to/user/HIPC/paper_developments/HIPC-Lung/exploratory_analyses/prol_assay_reps/jobs_scripts/jobs/run_prol_assay_reps_art-thold-0.5_2024-12-31/run_prol_assay_tcr_process.err.txt

module unload R
module load R/4.2.2

Rscript /path/to/user/HIPC/paper_developments/HIPC-Lung/exploratory_analyses/prol_assay_reps/jobs_scripts/prolif_assay_tcr_process.1.0.R --ArtifactThold 0.5
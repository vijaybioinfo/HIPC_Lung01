#!/bin/bash
#SBATCH --job-name=CompExtRef
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=80g
#SBATCH --time=160:00:00
#SBATCH --output=/path/to/user/HIPC/paper_developments/HIPC-Lung/ag-spc_tcr_reference/jobs_scripts/compile_external_ref_2024-06-19/run_compile_external_refs.out.txt
#SBATCH --error=/path/to/user/HIPC/paper_developments/HIPC-Lung/ag-spc_tcr_reference/jobs_scripts/compile_external_ref_2024-06-19/run_compile_external_refs.err.txt

module load R/3.6.1
Rscript /path/to/user/HIPC/paper_developments/HIPC-Lung/ag-spc_tcr_reference/jobs_scripts/compile_external_refs.1.0.R
#!/bin/bash
#SBATCH --job-name=ncrassa-nextflow
#SBATCH --partition=PARTITION_NAME
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=nextflow-%j.log

#Need to replace PARTITION_NAME above and below

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?This script must be submitted from the project directory}"

export NF_PARTITION="PARTITION_NAME"

# Need to replace below with cluster's module names

module load nextflow
module load sra-toolkit
module load fastqc

mkdir -p reports/execution

nextflow run main.nf \
	-profile slurm \
	--runs runs.txt \
	--outdir results \
	-work-dir PATH
	-with-report PATH
	-with trace PATH
	-resume

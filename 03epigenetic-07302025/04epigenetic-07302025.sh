#!/bin/bash
#SBATCH -p day # partition name
#SBATCH -t 10:00:00 # hours:minutes runlimit after which job will be killed
#SBATCH -c 1 # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 200G
#SBATCH --job-name Enmix_brain # Job name
#SBATCH -o %j.out # File to which standard out will be written
#SBATCH -e %j.err # File to which standard err will be written

# load conda
module load miniconda

# activate environment
conda activate ewas_ses_v04032024

# run Rscript for counting
Rscript 04epigenetic-07302025.Rscript

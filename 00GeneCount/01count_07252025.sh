#!/bin/bash
#SBATCH -p week # partition name
#SBATCH -t 2-03:00:00 # hours:minutes runlimit after which job will be killed
#SBATCH -c 1 # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 98G
#SBATCH --job-name Rsubread_counts # Job name
#SBATCH -o %j.out # File to which standard out will be written
#SBATCH -e %j.err # File to which standard err will be written

# load conda and activate environment
module load miniconda;conda activate rsubread_v04162024

# run Rscript for counting
Rscript 01count_07252025.Rscript 

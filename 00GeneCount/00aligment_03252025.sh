#!/bin/bash

#SBATCH -p week 		# partition name
#SBATCH -t 4-00:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -c 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 98G
#SBATCH --job-name STAR_aligment 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written

# Define paths
GENOME_DIR="/gpfs/gibbs/data/genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex/"
INPUT_DIRS=(
    "/gpfs/gibbs/pi/montalvo-ortiz/UTHealth_geno_transc/Transcriptomic/RNA2018"
    "/gpfs/gibbs/pi/montalvo-ortiz/UTHealth_geno_transc/Transcriptomic/RNA2023"
)
FASTQ_DIR="/vast/palmer/scratch/montalvo-ortiz/jjm262/00tclocks/00uthealth/00databases/00fastqs/"
STAR_OUTPUT_DIR="/vast/palmer/scratch/montalvo-ortiz/jjm262/00tclocks/00uthealth/00databases/01bams/"

# Create output directories if they don't exist
mkdir -p $FASTQ_DIR
mkdir -p $STAR_OUTPUT_DIR

# Number of threads
THREADS=6

# Step 1: Unzip all FASTQ files into the FASTQ_DIR
echo "Unzipping FASTQ files..."
for DIR in "${INPUT_DIRS[@]}"; do
    for FILE in $DIR/*.fq.gz; do
        BASENAME=$(basename "$FILE" .gz)  # Remove .gz
        DEST_FILE="${FASTQ_DIR}${BASENAME}"
        
        if [[ ! -f "$DEST_FILE" ]]; then
            echo "Unzipping $FILE to $FASTQ_DIR..."
            gunzip -c "$FILE" > "$DEST_FILE"
        else
            echo "File $DEST_FILE already exists, skipping..."
        fi
    done
done

# Step 2: Run STAR alignment using unzipped FASTQs
echo "Starting STAR alignment..."
for FILE in $FASTQ_DIR/*_1.fq; do
    SAMPLE=$(basename "$FILE" | sed 's/_1.fq//')

    READ1="${FASTQ_DIR}${SAMPLE}_1.fq"
    READ2="${FASTQ_DIR}${SAMPLE}_2.fq"

    # Run STAR
    STAR --genomeDir $GENOME_DIR \
         --runThreadN $THREADS \
         --readFilesIn $READ1 $READ2 \
         --outFileNamePrefix ${STAR_OUTPUT_DIR}${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         > ${STAR_OUTPUT_DIR}${SAMPLE}_STAR.log 2>&1

    echo "Finished processing $SAMPLE"
done

echo "Pipeline completed!"

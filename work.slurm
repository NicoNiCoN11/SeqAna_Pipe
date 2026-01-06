#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --output=logs/rna_seq_%j.out
#SBATCH --error=logs/rna_seq_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=24:00:00


# Create log directory
mkdir -p logs

# Load conda
source ~/.bashrc  # or wherever conda is initialized
conda activate rnaseq-pipeline

# Set work directory to scratch or large storage
export NXF_WORK="/home/jiguo/data/data/reverse_stranded/work"

# mkdir -p $NXF_WORK

# Run pipeline
nextflow run main.nf \
  --samplesheet /home/jiguo/SeqAna_Pipe/samplesheet.csv \
  --genome /home/jiguo/data/data/ref_genome/GRCh38.primary_assembly.genome.fa \
  --gtf /home/jiguo/data/data/ref_genome/gencode.v47.primary_assembly.annotation.gtf \
  --star_index /home/jiguo/data/data/outdir_rnaseq/star_index/star_index \
  --outdir /home/jiguo/data/data/reverse_stranded\
  -resume

# Cleanup work directory after successful completion (optional)
# Uncomment if you want automatic cleanup:
# if [ $? -eq 0 ]; then
#     echo "Pipeline successful, cleaning work directory..."
#     rm -rf $NXF_WORK
# fi
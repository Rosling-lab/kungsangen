#!/bin/bash

# runs pacbio RSII demultiplexing and base calling on a
# single node

#SBATCH -A snic2020-5-142
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 0-04:00:00
#SBATCH -J kungsängen
#SBATCH -C usage_mail
#SBATCH -M rackham
#SBATCH --mail-type=ALL
#SBATCH --output="logs/snakemake-%j.log"
#SBATCH --error="logs/snakemake-%j.log"

module load bioinfo-tools &&
module load snakemake &&

snakemake -pr --jobs $SLURM_JOB_CPUS_PER_NODE\
  --use-envmodules\
  --shadow-prefix /scratch\
  process/pb_363.ccs.bam\
  process/pb_363.laa.fastq

# Snakemake workflow file
# This workflow identifies raw sequencing reads and performs command-line based
# operations like circular consensus calling and demultiplexing.

import pandas as pd
import os.path
from glob import glob
import re
import subprocess
from snakemake.utils import listfiles

# For testing, parse the yaml file (this is automatically done by Snakemake)
#import yaml
#with open("config/config.yaml", 'r') as ymlfile: config = yaml.safe_load(ymlfile)

configfile: "config/config.yaml"

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

# load the dataset and region definitions
#datasets = pd.read_csv(config['dataset']).set_index('seq_run', drop = False)
#regions = pd.read_csv(config['regions']).set_index('region')

#### PacBio conversions ####
# find the PacBio movie files
moviefiles = [re.sub(r"\.bas\.h5", "", os.path.basename(m))
              for m in glob("raw/**/*.bas.h5", recursive = True)]

# endpoint target: convert all pacbio movies to Sequel format
localrules: convertmovies
rule convertmovies:
    input: expand("process/{movie}.{type}.bam",
               movie = moviefiles,
               type = ['subreads', 'scraps'])

# convert a raw RSII-format (.h5) movie to the Sequel format (.bam)
# these files are pretty large, so they are marked as temporary.
rule bax2bam:
    output:
        temp("process/{movie}.subreads.bam"),
        temp("process/{movie}.subreads.bam.pbi"),
        temp("process/{movie}.scraps.bam"),
        temp("process/{movie}.scraps.bam.pbi")
    input:
        lambda wildcards: glob("raw/**/{wildcards.movie}.*.bax.h5"
                               .format(wildcards = wildcards),
                               recursive = True)
    shadow: "shallow"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    group: "pacbio"
    params:
        prefix="process/{movie}"
    resources:
        walltime=10
    log: "logs/bax2bam_{movie}.log"
    shell: "bax2bam {input} -o {params.prefix} &> {log}"

# merge all the movies which belong to the same plate
rule mergebam:
    output: temp("process/pb_363.subreads.bam")
    input:
        expand("process/{movie}.subreads.bam",
               movie = moviefiles[wildcards.seqplate])
    shadow: "shallow"
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    group: "pacbio"
    resources:
        walltime=5
    log: "logs/mergebam.log"
    shell: "samtools merge {output} {input} &>{log}"

# demultiplex pacbio subreads using lima
rule lima:
    output:
        temp("process/pb_363.demux.subreads.bam")
    input:
        bam = "process/pb_363.subreads.bam",
        tags = "tags/its1_lr5_barcodes.fasta"
    shadow: "shallow"
    conda: "conda/pacbiodemux.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    group: "pacbio"
    resources:
        walltime=20
    log: "logs/lima.log"
    shell: "lima {input.bam} {input.tags} {output} --different --peek-guess &>{log}"

# find haplotypes (ASVs) from pacbio subreads
rule laa:
    output:
        result="process/pb_363.laa.fastq",
        junk="process/pb_363.chimera_noise.fastq",
        report="process/pb_363.report.csv",
        pcr="process/pb_363.pcr.csv"
    input: "process/pb_363.demux.subreads.bam"
    shadow: "shallow"
    conda: "conda/pacbiolaa.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    group: "pacbio"
    params:
        prefix="process/pb_363"
    resources:
        walltime=240
    log: "logs/laa_{{seqplate}}.log".format_map(config)
    shell:
        """
        laa {input}\\
            --maxReads 10000\\
            --maxClusteringReads 10000\\
            --minLength 1000\\
            --minClusterSize 3\\
            --resultFile {output.result}\\
            --junkFile {output.junk}\\
            --reportFile {output.report}\\
            --inputReportFile {output.pcr}\\
            --subreadsReportPrefix {output.prefix}
        """
        

# generate a circular consensus sequence from raw PacBio reads
rule ccs:
    output: temp("process/pb_363.ccs.bam")
    input: "process/pb_363.demux.subreads.bam"
    resources:
        walltime=120
    shadow: "shallow"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    group: "pacbio"
    threads: 4
    log: "logs/ccs.log"
    shell:
         "ccs --numThreads {threads} {input} {output} &>{log}"

# Snakemake workflow file
# This workflow identifies raw sequencing reads and performs command-line based
# operations like circular consensus calling and demultiplexing.

import os.path
from glob import glob
import re
import subprocess
from math import gcd

# For testing, parse the yaml file (this is automatically done by Snakemake)
#import yaml
#with open("config/config.yaml", 'r') as ymlfile: config = yaml.safe_load(ymlfile)

#configfile: "config/config.yaml"

# Find the maximum number of cores available to a single node on SLURM,
# or if we aren't in a SLURM environment, how many we have on the local machine.
try:
    maxthreads = max([int(x) for x in re.findall(r'\d+', subprocess.check_output(["sinfo", "-O", "cpus"]).decode())])
except FileNotFoundError:
    maxthreads = int(subprocess.check_output("nproc").decode())

# find the PacBio movie files
moviefiles = [re.sub(r"\.bas\.h5", "", os.path.basename(m))
              for m in glob("raw/**/*.bas.h5", recursive = True)]

# figure out how to distribute the movie files equally over the cores we have available.
moviejobs = gcd(maxthreads, len(moviefiles))
moviethreads = maxthreads/moviejobs

localrules: all
rule all:
    input:
        "process/pb_363.ccs.fastq.gz",
        "process/pb_363.laa.fastq.gz",
        "process/pb_363.laagc.fastq.gz",
        "process/pb_363.ccs.swarm.cons.fasta"

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
    params:
        prefix="process/{movie}"
    resources:
        walltime=20
    threads: 2
    log: "logs/bax2bam_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    shell: "bax2bam {input} -o {params.prefix} &> {log}"

# endpoint target: convert all pacbio movies to Sequel format
rule convertmovies:
    input:
        expand("process/{movie}.{type}.bam",
               movie = moviefiles,
               type = ['subreads', 'scraps'])

# merge all the movies which belong to the same plate
rule mergebam:
    output: temp("process/pb_363.subreads.bam")
    input:
        expand("process/{movie}.subreads.bam",
               movie = moviefiles)
    shadow: "shallow"
    resources:
        walltime=20
    log: "logs/mergebam.log"
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    shell: "samtools merge {output} {input} &>{log}"

# a seqrun identifier should start with pb (RSII) or ps (Sequel), then 3 digits,
# and optionally 3 more if multiple plates were sequenced
wildcard_constraints:
    seqrun = "p[bs]_\d{3}(_\d{3})?",
    movie = "m\d{6}_.*"

# demultiplex pacbio subreads using lima
# the output is still a single BAM file, but it has the barcode assignments
# in the headers for each sequence
rule lima:
    output:
        temp("process/{movie}.subreads.demux.bam")
    input:
        bam = "process/{movie}.subreads.bam",
        tags = "tags/its1_lr5_barcodes.fasta"
    shadow: "shallow"
    threads: moviethreads
    resources:
        walltime=20
    log: "logs/lima_{movie}.log"
    conda: "conda/pacbiodemux.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "lima {input.bam} {input.tags} --different --peek-guess -j {threads} --keep-tag-idx-order {output} &>{log}"

# filter out the samples which are not being used in this project.
rule sieve:
    output:
          "process/{movie}.subreads.demux.sieve.bam"
    input:
         bam="process/{movie}.subreads.demux.bam",
         samples = "tags/which_tags.txt"
    shadow: "shallow"
    threads: 1
    resources:
        walltime=5
    log: "logs/sieve_{movie}.log"
    conda: "conda/pacbiodemux.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "bamsieve --barcodes --whitelist {input.samples} {input.bam} {output} &>{log}"

# endpoint target: demultiplex and sieve all movies
rule sievemovies:
    input:
        expand("process/{movie}.subreads.demux.sieve.bam",
               movie = moviefiles)

# generate circular consensus sequences from subreads
rule ccs:
    output: "process/{movie}.ccs.bam"
    input: "process/{movie}.subreads.demux.sieve.bam"
    resources:
        walltime=120
    shadow: "shallow"
    threads: moviethreads
    log: "logs/ccs_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    shell: "ccs --numThreads {threads} {input} {output} &>{log}"

# convert a ccs BAM to a fastq
# this loses a lot of PacBio-specific information, but it is useful for other software.
rule bam2fastq:
    output: temp("process/{movie}.ccs.fastq.gz")
    input: "process/{movie}.ccs.bam"
    resources:
             walltime=10
    threads: 1
    log: "logs/bam2fastq_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "smrt/5.0.1"
    shell: "bam2fastq -o process/{wildcards.movie}.ccs {input} &>{log}"

# quality filter the ccs and dereplicate
# allow up to 15 expected errors (about 1%) and minimum length 1000
# the output has only one entry for each unique sequence
# the label is the ZMW name of one of the first appearance, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output:
        fasta="process/pb_363.ccs.derep.fasta",
        uc="process/pb_363.ccs.derep.uc"
    input: expand("process/{movie}.ccs.fastq.gz", movie = moviefiles)
    resources:
        walltime=10
    shadow: "shallow"
    threads: 2
    log: "logs/derep_pb_363.log"
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1"
    shell:
        """
         zcat {input} |
         vsearch --fastq_filter - \\
            --fastq_maxee 15 \\
            --fastq_qmax 93 \\
            --fastq_minlen 1000 \\
            --fastaout - |
         vsearch --derep_fulllength - \\
            --sizeout \\
            --fasta_width 0\\
            --output {output.fasta}\\
            --uc {output.uc}    
        """

# Swarm-cluster the CCS reads
# this is basically the same thing as single-linkage clustering
# we use a maximum distance of 30, which is double the max ee
# so theoretically, every sequence which came from a single biological variant
# should end up in the same cluster.
rule gefast:
    output: "process/{seqrun}.ccs.swarm"
    input: "process/{seqrun}.ccs.derep.fasta"
    resources:
             walltime=120
    shadow: "shallow"
    group: "pacbio"
    threads: 8
    log: "logs/gefast_{seqrun}.log"
    conda: "conda/gefast.yaml"
    # not available as a module on UPPMAX
    shell:
         """
         GeFaST {input} \\
            --threshold 30\\
            --sep-abundance ";size="\\
            --swarm-output {output}\\
            --swarm-fastidious\\
            --swarm-no-otu-breaking\\
            --swarm-num-explorers {threads}\\
            --swarm-num-grafters {threads}\\
            --swarm-num-threads-per-check {threads} &>{log}
         """

rule swarmselect:
    output: directory("process/swarm/{seqrun}")
    input:
        swarm="process/{seqrun}.ccs.swarm",
        uc=   "process/{seqrun}.ccs.derep.uc",
        bam=  "process/{seqrun}.ccs.bam"
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    shell:
        """
        cat {input.swarm} |
          parallel --pipe -N1 \\
          "tr 
        """


# split the swarm file up
rule swarmconsensus:
    output:
        consensus="process/{seqrun}.ccs.swarm.cons.fasta",
        swarmdir=directory("process/swarm/{seqrun}")
    input:
        swarm="process/{seqrun}.ccs.swarm",
        script="scripts/swarm_consensus.sh",
        fasta="process/{seqrun}.ccs.derep.fasta"
    resources:
             walltime=240
    threads: maxthreads
    conda: "conda/cons.yaml"
    shell:
        """
        [ -d {output.swarmdir} ] || mkdir -p {output.swarmdir}
        cat {input.swarm} |
        parallel --pipe -j{threads} -N1 {input.script} {input.fasta} {output.swarmdir} {{#}} {wildcards.seqrun}
        ls {output.swarmdir}/swarm_*.cons.fasta | xargs cat >{output.consensus}
        """

# find haplotypes (ASVs) from pacbio subreads using gefast cluster consensus as guides
rule laagc:
    output:
        result="process/{seqrun}.laagc.fastq.gz",
        junk="process/{seqrun}.chimera_noisegc.fastq.gz",
        report="process/{seqrun}.reportgc.csv",
        pcr="process/{seqrun}.pcrgc.csv"
    input:
        subreads="process/{seqrun}.demux.subreads.bam",
        guide="process/{seqrun}.ccs.swarm.cons.fasta"
    shadow: "shallow"
    threads: maxthreads
    params:
        prefix="process/pb_363",
        result_prefix = "process/pb_363.laa.fastq",
        junk_prefix = "process/pb_363.chimera_noise.fastq"
    resources:
        walltime=240
    log: "logs/laa_{seqrun}.log"
    conda: "conda/pacbiolaa.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    shell:
        """
        laa {input.subreads}\\
            --numThreads {threads}\\
            --clusterGuide {input.guide}\\
            --ignoreBc\\
            --maxReads 100000\\
            --maxClusteringReads 20000\\
            --minLength 1000\\
            --minClusterSize 3\\
            --maxPhasingReads 20000\\
            --minSplitReads 3\\
            --minSplitFraction 0.01\\
            --resultFile {params.result_prefix}\\
            --junkFile {params.junk_prefix}\\
            --reportFile {output.report}\\
            --inputReportFile {output.pcr}\\
            --subreadsReportPrefix {params.prefix} >&{log} && 
        gzip {params.result_prefix} >&{log} &&
        gzip {params.junk_prefix} >&{log}
        """

# find haplotypes (ASVs) from pacbio subreads using gefast cluster consensus as guides
rule laa:
    output:
        result="process/{seqrun}.laa.fastq.gz",
        junk="process/{seqrun}.chimera_noise.fastq.gz",
        report="process/{seqrun}.report.csv",
        pcr="process/{seqrun}.pcr.csv"
    input:
        subreads="process/{seqrun}.demux.subreads.bam"
    shadow: "shallow"
    threads: maxthreads
    params:
        prefix="process/pb_363",
        result_prefix = "process/pb_363.laa.fastq",
        junk_prefix = "process/pb_363.chimera_noise.fastq"
    resources:
        walltime=240
    log: "logs/laa_{seqrun}.log"
    conda: "conda/pacbiolaa.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    shell:
        """
        laa {input.subreads}\\
            --numThreads {threads}\\
            --ignoreBc\\
            --maxReads 100000\\
            --maxClusteringReads 20000\\
            --minLength 1000\\
            --minClusterSize 3\\
            --maxPhasingReads 20000\\
            --minSplitReads 3\\
            --minSplitFraction 0.01\\
            --resultFile {params.result_prefix}\\
            --junkFile {params.junk_prefix}\\
            --reportFile {output.report}\\
            --inputReportFile {output.pcr}\\
            --subreadsReportPrefix {params.prefix} >&{log} && 
        gzip {params.result_prefix} >&{log} &&
        gzip {params.junk_prefix} >&{log}
        """
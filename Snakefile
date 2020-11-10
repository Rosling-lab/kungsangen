# Snakemake workflow file
# This workflow identifies raw sequencing reads and performs command-line based
# operations like circular consensus calling and demultiplexing.

import os.path
from glob import glob
import re
import subprocess

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

# load the dataset and region definitions
#datasets = pd.read_csv(config['dataset']).set_index('seq_run', drop = False)
#regions = pd.read_csv(config['regions']).set_index('region')

#### PacBio conversions ####
# find the PacBio movie files
moviefiles = [re.sub(r"\.bas\.h5", "", os.path.basename(m))
              for m in glob("raw/**/*.bas.h5", recursive = True)]

localrules: all
rule all:
    input:
        "process/pb_363.ccs.fastq.gz",
        "process/pb_363.laa.fastq.gz",
        "process/pb_363.laagc.fastq.gz",
        "process/pb_363.ccs.swarm.cons.fasta"

# endpoint target: convert all pacbio movies to Sequel format
rule convertmovies:
    input:
        expand("process/{movie}.{type}.bam",
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
        temp("process/{seqrun}.demux.subreads.bam")
    input:
        bam = "process/{seqrun}.subreads.bam",
        tags = "tags/its1_lr5_barcodes.fasta"
    shadow: "shallow"
    threads: maxthreads
    resources:
        walltime=20
    log: "logs/lima_{seqrun}.log"
    conda: "conda/pacbiodemux.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "lima {input.bam} {input.tags} {output} --different --peek-guess &>{log} -j {threads}"

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

# generate circular consensus sequences from raw PacBio reads
# this preserves the demultiplexing info in the headers
rule ccs:
    output: "process/{seqrun}.ccs.bam"
    input: "process/{seqrun}.demux.subreads.bam"
    resources:
        walltime=120
    shadow: "shallow"
    threads: maxthreads
    log: "logs/ccs_{seqrun}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1"
    shell: "ccs --numThreads {threads} {input} {output} &>{log}"

# convert a ccs BAM to a fastq
# this loses a lot of PacBio-specific information, but it is useful for other software.
rule bam2fastq:
    output: "process/{seqrun}.ccs.fastq.gz"
    input: "process/{seqrun}.ccs.bam"
    resources:
             walltime=10
    shadow: "shallow"
    threads: 1
    log: "logs/bam2fastq_{seqrun}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "smrt/5.0.1"
    shell: "bam2fastq -o process/{wildcards.seqrun}.ccs -p ccs {input} &>{log}"

# quality filter the fastq and dereplicate it
# allow up to 15 expected errors (about 1%) and minimum length 1000
# the output has only one entry for each unique sequence
# the label is the md5 hash of the sequence, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output: "process/{seqrun}.ccs.derep.fasta"
    input: "process/{seqrun}.ccs.fastq.gz"
    resources:
        walltime=10
    shadow: "shallow"
    threads: 2
    log: "logs/derep_{seqrun}.log"
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1"
    shell:
        """
         vsearch --fastq_filter {input}\\
            --fastq_maxee 15 \\
            --fastq_qmax 93 \\
            --fastq_minlen 1000 \\
            --fastaout - |
         vsearch --derep_fulllength - \\
            --relabel_md5 \\
            --sizeout \\
            --fasta_width 0\\
            --output {output}    
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
            --swarm-num-threads-per-check {threads}
         """

# split the swarm file up into a separate file for each cluster
# these files give the headers of the sequences which belong to each cluster
checkpoint splitswarm:
    output: directory("process/swarm/{seqrun}")
    input: "process/{seqrun}.ccs.swarm"
    resources:
             walltime=5
    threads: 1
    log: "logs/splitswarm_{seqrun}.log"
    shell:
        """
        [ -d {output} ] || mkdir -p {output} &&
        split -l 1 -a 5 --numeric-suffixes {input} {output}/swarm_ &>{log}
        """

# extract the sequences which belong to a cluster from the dereplicated fasta file
rule extractswarm:
    output: "process/swarm/{seqrun}/swarm_{num}.fasta"
    input:
        swarm="process/swarm/{seqrun}/swarm_{num}",
        fasta="process/{seqrun}.ccs.derep.fasta"
    resources:
        walltime=1
    threads: 1
    log: "logs/extractswarm_{seqrun}_{num}.log"
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1"
    shell:
        """
        tr " " "\\n" {input.swarm} |
        vsearch --fastx_getseqs {input.fasta} \\
            --labels -\\
            --fastaout {output}
        
        """

# align the sequences in each swarm cluster
# this makes a fast alignment in MUSCLE, because the sequences are all VERY similar
# The gap penalty is small, because PacBio has lots of small indel errors.
rule align:
    output: "process/swarm/{seqrun}/swarm_{num}.aln.fasta"
    input: "process/swarm/{seqrun}/swarm_{num}.fasta"
    resources:
        walltime=30
    threads: 1
    log: "logs/alignswarm_{seqrun}_{num}.log"
    conda: "conda/muscle.yaml"
    envmodules:
        "bioinfo-tools",
        "muscle/3.8.1551"
    shell: "muscle -in {input} -out {output} -maxiters 1 -diags -gapopen -0.5"

# Calculate the consensus sequence for each swarm cluster
# First we need to re-replicate the file
# tools like VSEARCH and fastx_toolkit won't accept sequences with gaps, so we
# do this with sed, tr, and awk
rule consensus:
    output: "process/swarm/{seqrun}/swarm_{num}.cons.fasta"
    input: "process/swarm/{seqrun}/swarm_{num}.aln.fasta"
    threads: 1
    log: "logs/consensus_{seqrun}_{num}.log"
    conda: "conda/cons.yaml"
    envmodules:
        "bioinfo-tools",
        "emboss/6.6.0"
    shell:
         """         
         # linearize the fasta; i.e. put each sequence on one line, with tab
         # between header and seq.
         sed '/^>/s/$/@/; s/^>/#>/' {input} |
         tr -d "\n" |
         tr "#" "\n" |
         tr "@" "\t" |
         # duplicate each line based on ";size=" from the header
         gawk '{c=gensub(/.+;size=([0-9]+).*/, "\\1", 1, $1); c=int(c); while (c--) {print $1; print $2;}}' |
         # calculate the consensus
         cons -filter -sformat fasta -osformat fasta -plurality 0.3 -name swarm{wildcards.num} |
         sed '/^>/!s/n//g' >{output} 2>{log}
         """

#
def gather_swarm(wildcards):
    checkpoint_output = checkpoints.splitswarm.get(**wildcards).output[0]
    return expand("process/swarm/{seqrun}/swarm_{i}.cons.fasta",
                  seqrun = wildcards.seqrun,
                  i=glob_wildcards(os.path.join(checkpoint_output, "swarm_{i,\d+}")).i)

# Gather the consensuses from all the clusters into one fasta file
rule gather:
    output: "process/{seqrun}.ccs.swarm.cons.fasta"
    input: gather_swarm
    shell: "cat {input} >{output}"
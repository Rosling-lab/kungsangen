# Snakemake workflow file
# This workflow identifies raw sequencing reads and performs command-line based
# operations like circular consensus calling and demultiplexing.

import os.path
from glob import glob
import re
import subprocess
from math import gcd
from snakemake.io import glob_wildcards

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
        "SMRT/5.0.1" # no bax2bam in newer versions
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
        tags = "tags/fwd_rev_barcodes.fasta"
    shadow: "shallow"
    threads: moviethreads
    resources:
        walltime=20
    log: "logs/lima_{movie}.log"
    conda: "conda/pacbiodemux.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "lima {input.bam} {input.tags} --different --peek-guess -j {threads} {output} &>{log}"

# filter out the samples which are not being used in this project.
# bamsieve from pacbio looked like it would be a good way to do this,
# but it doesn't let you supply two barcodes for asymmetric barcoding
# samtools view doesn't filter by the bc: tag
# so grep to the rescue!
rule sieve:
    output:
          "process/{movie}.subreads.demux.sieve.bam"
    input:
         bam="process/{movie}.subreads.demux.bam",
         samples = "tags/which_tags.txt"
    params:
        samples = "{movie}.tags.txt"
    shadow: "shallow"
    threads: 4
    resources:
        walltime=20
    log: "logs/sieve_{movie}.log"
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    shell:
        """
        sed 's/^/bc:B:S,/' {input.samples} > {params.samples} 2>{log}
        {{
            samtools view -H {input.bam}
            samtools view -@1 {input.bam} |
            grep -f {params.samples}
        }} |
        samtools view -@1 - -o {output} 2>>{log}
        rm {params}
        """

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
        "SMRT/5.0.1" # ccs from newer versions doesn't accept RSII data
    shell: "ccs --numThreads {threads} --richQVs {input} {output} &>{log}"

rule ccs2:
    output: "process/nodemux/{movie}.ccs.bam"
    input: "process/{movie}.subreads.bam"
    resources:
        walltime=120
    shadow: "shallow"
    threads: moviethreads
    log: "logs/nodemux_ccs_{movie}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" # ccs from newer versions doesn't accept RSII data
    shell: "[ -d process/nodemux ] || mkdir -p process/nodemux && ccs --numThreads {threads} {input} {output} &>{log}"


# convert a ccs BAM to a fastq
# this loses a lot of PacBio-specific information, but it is useful for other software.
rule bam2fastq:
    output: temp("process/{name}.fastq.gz")
    input:
        bam = "process/{name}.bam",
        pbi = "process/{name}.bam.pbi"
    resources:
             walltime=10
    threads: 1
    log: "logs/bam2fastq_{name}.log"
    conda: "conda/pacbio.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "bam2fastq -o process/{wildcards.name} {input.bam} &>{log}"

rule nodemux_combine:
    output: "process/nodemux/pb_363.ccs.fastq.gz"
    input: expand("process/nodemux/{movie}.ccs.fastq.gz", movie = moviefiles)
    resources:
             walltime=10
    threads: 1
    log: "logs/nodemux_ccs_pb_363.log"
    shell: "cat {input} >{output} 2>{log}"

# lima doesn't store any information about orientation when run on subreads,
# and the primers are already gone, so we can't orient using primers.
# instead, search for 5.8S.
rule orient:
    output:
        orient = "process/{movie}.ccs.orient.fastq.gz",
        noprimer = "process/{movie}.ccs.noprimer.fastq.gz"
    input:
        ccs = "process/{movie}.ccs.fastq.gz"
    threads: moviethreads
    log: "logs/{movie}_orient.log"
    conda: "conda/orient.yaml"
    shell:
        """
        cutadapt\\
            -a "TCCGTAGGTGAACCTGC;e=0.15...CGAAGTTTCCCTCAGGA;required;e=0.15"\\
            --action=none\\
            --revcomp\\
            -o {output.orient}\\
            --untrimmed-output {output.noprimer}\\
            -j {threads}\\
            {input.ccs}
        """

# quality filter the ccs and dereplicate
# allow up to 15 expected errors (about 1%) and minimum length 1000
# the output has only one entry for each unique sequence
# the label is the ZMW name of one of the first appearance, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output:
        fasta="process/pb_363.ccs.derep.fasta",
        fastq="process/pb_363.ccs.orient.fastq.gz",
        uc="process/pb_363.ccs.derep.uc"
    input: expand("process/{movie}.ccs.orient.fastq.gz", movie = moviefiles)
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
         fastq=$(mktemp --suffix .fastq) &&
         trap 'rm ${{fastq}}' EXIT &&
         zcat {input} |
         vsearch --fastq_filter - \\
            --fastq_maxee 15 \\
            --fastq_qmax 93 \\
            --fastq_minlen 1000 \\
            --fastqout ${{fastq}}\\
            --fastaout - |
         vsearch --derep_fulllength - \\
            --sizeout \\
            --fasta_width 0\\
            --output {output.fasta}\\
            --uc {output.uc} &&
         gzip -c ${{fastq}} >{output.fastq}
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

# index a pacbio bam file
rule index:
    output: "{basename}.bam.pbi"
    input: "{basename}.bam"
    conda: "conda/swarmextract.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "pbindex {input}"

# for each of the swarm clusters, make a BAM file containing the source CCS
checkpoint swarmselect:
    output: directory("process/swarm/{seqrun}_{type}")
    input:
        swarm="process/{seqrun}.ccs.swarm",
        uc=   "process/{seqrun}.ccs.derep.uc",
        bam=  expand("process/{movie}.{{type}}.bam", movie = moviefiles),
        pbi= expand("process/{movie}.{{type}}.bam.pbi", movie = moviefiles),
        script="scripts/swarm_laa.sh"
    log: "logs/swarmselect_{seqrun}_{type}.log"
    threads: maxthreads
    conda: "conda/swarmextract.yaml"
    #bamsieve in the env module gives mysterious errors
    #envmodules:
    #    "bioinfo-tools",
    #    "samtools",
    #    "SMRT/7.0.1",
    #    "gnuparallel/20180822"
    shell:
        """
	[ -d {output} ] || mkdir -p {output}
        cat {input.swarm} |
          {{ parallel --pipe -N1 -j {threads} {input.script} {{#}} {input.uc} process .{wildcards.type} {output}; }} &>{log}
        """

# get consensus of CCS reads from each swarm cluster
rule swarm_consensus:
    output: "process/{seqrun}.swarm.cons.fasta"
    input:
        swarm="process/{seqrun}.ccs.swarm",
        uc=   "process/{seqrun}.ccs.derep.uc",
        fastq=  "process/{seqrun}.ccs.orient.fastq.gz",
        script="scripts/swarm_consensus.sh",
        c3s= "bin/c3s"
    log: "logs/swarm_consensus_{seqrun}.log"
    threads: maxthreads
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1",
        "gnuparallel/20180822"
    shell:
        """
        cat {input.swarm} |
        {{ parallel --pipe -N1 -j {threads} {input.script} {wildcards.seqrun} {{#}} {input.uc} {input.fastq}; }} >{output}
        """

rule table_translate:
    output: "process/sample.tsv"
    input: "start_files/new_fqNames.txt"
    log: "logs/table_translate.log"
    threads: 1
    shell:
        """
        awk -F '[-,_.]' '{{print $2-1 "," $4-1, $6  "_"  $7  "_"  $8}}' {input} >{output} 2>{log}
        """

rule swarm_table:
    output: "process/{seqrun}_{type}.swarm.table"
    input:
        bamdir = "process/swarm/{seqrun}_{type}",
        script = "scripts/swarm_table.sh",
        samples = "process/sample.tsv"
    log: "logs/swarm_table_{seqrun}_{type}.log"
    threads: maxthreads
    conda: "conda/samtools_parallel.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools/1.10",
        "gnuparallel/20180822"
    shell:
        """
        ls -1 {input.bamdir}/swarm_*.bam | {{ parallel -j {threads} {input.script} {{}} {input.samples}; }} >{output} 2>{log}
        """

# find haplotypes (ASVs) from pacbio subreads or ccs in each gefast cluster
rule laa:
    output:
        result="process/swarm/{seqrun}_{type}/swarm_{clust}.laa.fastq.gz",
        junk="process/swarm/{seqrun}_{type}/swarm_{clust}.junk.fastq.gz",
        report="process/swarm/{seqrun}_{type}/swarm_{clust}.report.csv",
        subreads="process/swarm/{seqrun}_{type}/swarm_{clust}.subreads.csv",
        pcr="process/swarm/{seqrun}_{type}/swarm_{clust}.pcr.csv"
    input:
        subreads="process/swarm/{seqrun}_{type}/swarm_{clust}.bam"
    shadow: "shallow"
    # although each invocation of laa can run in parallel, this is probably quite inefficient for the
    # majority of clusters, which are small.  So just do single threads, with multiple clusters running in parallel.
    threads: 1
    params:
        clustersize=lambda wildcards : '3' if wildcards.type == "ccs" else '50'
    resources:
        walltime=240
    log: "logs/laa_{seqrun}_{type}_{clust}.log"
    conda: "conda/pacbiolaa.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/5.0.1" #laa from newer versions doesn't accept RSII data
    shell:
        """
        laa {input.subreads}\\
            --numThreads {threads}\\
            --ignoreBc\\
            --maxReads 100000\\
            --maxClusteringReads 10000\\
            --minLength 1000\\
            --maxLength 2500\\
            --minClusterSize 3\\
            --maxPhasingReads 100000\\
            --minSplitReads 3\\
            --minSplitFraction 0.001\\ >&{log} && 
        touch amplicon_analysis.fastq &&
        touch amplicon_analysis_chimeras_noise.fastq &&
        touch amplicon_analysis_summary.csv &&
        touch amplicon_analysis_input.csv &&
        touch amplicon_analysis_subreads.csv &&
        gzip -c amplicon_analysis.fastq >{output.result} 2>>{log} &&
        gzip -c amplicon_analysis_chimeras_noise.fastq >{output.junk} 2>>{log} &&
        mv amplicon_analysis_summary.csv {output.report} &>>{log} &&
        mv amplicon_analysis_input.csv {output.pcr} &>>{log} &&
        mv amplicon_analysis_subreads.csv {output.subreads} &>>{log}
        """

def swarmfiles(wildcards):
    swarmdir = checkpoints.swarmselect.get(**wildcards).output[0]
    clusters = glob_wildcards(os.path.join(swarmdir, "swarm_{cluster}.bam"))
    return expand("{swarmdir}/swarm_{cluster}.{ext}", swarmdir = swarmdir, cluster = clusters.cluster,
                  ext = ["laa.fastq.gz", "junk.fastq.gz", "report.csv", "subreads.csv", "pcr.csv"])
rule all_laa:
    output: touch("process/all_laa_{seqrun}_{type}")
    input: swarmfiles

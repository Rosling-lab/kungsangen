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
        #"process/pb_363.laa.fastq.gz",
        #"process/pb_363.laagc.fastq.gz",
        "process/pb_363.swarm.cons.fasta",
        "process/pb_363.vclust.cons.fasta",
        "process/pb_363_ccs.swarm.table",
        "process/pb_363_ccs.vclust.table"

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
        expand("process/{{movie}}.subreads.demux.lima.{ext}", ext = ["clips", "counts", "guess", "report", "summary"]),
        bam = temp("process/{movie}.subreads.demux.bam")
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
    shell: "lima {input.bam} {input.tags} --different --peek-guess -j {threads} {output.bam} &>{log}"

# sort the lima report alphabetically
# by default it is sorted numerically
rule sortreport:
    output: "process/{infile}.sort.lima.report"
    input: "process/{infile}.lima.report"
    shell: "tail -n +1 {input} | sort -k 1b,1 >{output}"

# count the number of reads, unique reads, and zmws in a BAM file
rule bamstats:
    output: "process/{base}.bamstats"
    input: "process/{base}.bam"
    wildcard_constraints:
        base = "m.+"
    threads: moviethreads
    params: samthreads = lambda wildcards, threads: threads - 1
    conda: "conda/samtools.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools"
    shell:
        """
        samtools view -@{params.samthreads} {input} |
        gawk -F'[/\\t ]' '{{lc++; count[$1,$2]++; derep[$12]++}}; END{{print lc, length(derep), length(count)}}' >{output}
        """

rule allbamstats:
    output: "process/all.bamstats"
    input:
        rawbamstats = expand("process/{movie}.subreads.bamstats", movie = moviefiles),
        sievebamstats = expand("process/{movie}.subreads.demux.sieve.bamstats", movie = moviefiles),
        ccsbamstats = expand("process/{movie}.ccs.bamstats", movie = moviefiles),
        limastats = expand("process/{movie}.subreads.demux.lima.report", movie = moviefiles)
    shell:
        """
        echo "step\treads\tzmw\tunique" >{output}
        gawk 'BEGIN{{OFS="\\t"}}; {{reads+=$1; zmw+=$3; unique+=$2}}; END{{print "raw", reads, zmw, unique}}' {input.rawbamstats} >>{output}
        gawk 'BEGIN{{OFS="\\t"}}; FNR>1{{reads+=$28+2; zmw++; unique=reads}}; END{{print "demux", reads, zmw, unique}}' {input.limastats} >>{output}
        gawk 'BEGIN{{OFS="\\t"}}; {{reads+=$1; zmw+=$3; unique+=$2}}; END{{print "sieve", reads, zmw, unique}}' {input.sievebamstats} >>{output}
        gawk 'BEGIN{{OFS="\\t"}}; {{reads+=$1; zmw+=$3; unique+=$2}}; END{{print "ccs", reads, zmw, unique}}' {input.ccsbamstats} >>{output}
        """

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
# so orient using the primers
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
            -a "TCCGTAGGTGAACCTGC;e=0.15;o=10...CGAAGTTTCCCTCAGGA;required;e=0.15;o=10"\\
            --action=none\\
            --revcomp\\
            -o {output.orient}\\
            --untrimmed-output {output.noprimer}\\
            -j {threads}\\
            {input.ccs}
        """

# quality filter the ccs and dereplicate
# allow up to 1% expected errors, minimum length 1000, maximum length 2000
# the fasta output has only one entry for each unique sequence
# the label is the ZMW name of one of the first appearance, followed by ";size=n"
# "n" gives the number of times the sequence appears.
rule derep:
    output:
        fasta="process/pb_363.ccs.derep.fasta",
        fastq="process/pb_363.ccs.orient.fastq.gz",
        trimmed="process/pb_363.ccs.trimmed.fastq.gz",
        tooshort="process/pb_363.ccs.orient.tooshort.fastq.gz",
        toolong="process/pb_363.ccs.orient.toolong.fastq.gz",
        toopoor="process/pb_363.ccs.orient.toopoor.fastq.gz",
        uc="process/pb_363.ccs.derep.uc"
    input: expand("process/{movie}.ccs.orient.fastq.gz", movie = moviefiles)
    resources:
        walltime=10
    shadow: "shallow"
    threads: 2
    log: "logs/derep_pb_363.log"
    conda: "conda/orient.yaml"
    shell:
        """
         fastq=$(mktemp --suffix .fastq) &&
         tooshort=$(mktemp --suffix .fastq) &&
         toolong=$(mktemp --suffix .fastq) &&
         toopoor=$(mktemp --suffix .fastq) &&
         trap 'rm ${{fastq}} ${{tooshort}} ${{toolong}} ${{toopoor}}' EXIT &&
         cat {input} > {output.fastq} &&
         cutadapt \\
            -a "TCCGTAGGTGAACCTGC;e=0.15;o=10...CGAAGTTTCCCTCAGGA;required;e=0.15;o=10"\\
            -o -\\
            -j 1\\
            {output.fastq} |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_maxee 12 \\
            --fastq_qmax 93 \\
            --fastqout_discarded ${{toopoor}} \\
            --fastqout - |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_minlen 50 \\
            --fastqout_discarded ${{tooshort}} \\
            --fastqout - |
         vsearch --fastq_filter - \\
            --threads 1 \\
            --fastq_qmax 93 \\
            --fastq_maxlen 2999 \\
            --fastqout_discarded ${{toolong}} \\
            --fastqout ${{fastq}}\\
            --fastaout - |
         vsearch --derep_fulllength - \\
            --threads 1 \\
            --sizeout \\
            --fasta_width 0\\
            --output {output.fasta} \\
            --uc {output.uc} &&
         gzip -c ${{fastq}} >{output.trimmed} &&
         gzip -c ${{tooshort}} >{output.tooshort} &&
         gzip -c ${{toolong}} >{output.toolong} &&
         gzip -c ${{toopoor}} >{output.toopoor}
        """

# Centroid-cluster the CCS reads
rule vclust:
    output:
        uc = "process/{seqrun}.ccs.vclust.uc",
        fasta = "process/{seqrun}.ccs.vclust.fasta",
        otutab = "process/{seqrun}.ccs.vclust.otu_table.txt",
        clustfile = "process/{seqrun}.ccs.vclust"
    input: "process/{seqrun}.ccs.derep.fasta"
    params:
        clusterdir = "process/vclust/{seqrun}"
    shadow: "shallow"
    threads: maxthreads
    log: "logs/vclust_{seqrun}.log"
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1"
    shell:
        """
        mkdir -p {params.clusterdir} &&
        vsearch --cluster_size {input} \\
            --sizein \\
            --consout {output.fasta} \\
            --uc {output.uc} \\
            --otutabout {output.otutab} \\
            --id 0.99 \\
            --clusterout_id \\
            --threads {threads} \\
            --clusters {params.clusterdir}/otu &&
        for clust in $(ls {params.clusterdir}); do
            sed -n '/^>/s/^>//p' <{params.clusterdir}/$clust | tr "\n" " " &&
            echo
        done |
        grep -v "^$" >{output.clustfile} &&
        rm -r {params.clusterdir}
        """


# Swarm-cluster the CCS reads
# this is basically the same thing as single-linkage clustering
# we use a maximum distance of 30, which is ~double the max ee
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
            --swarm-fastidious 0\\
            --swarm-no-otu-breaking\\
            --swarm-num-explorers {threads}\\
            --swarm-num-grafters {threads}\\
            --swarm-num-threads-per-check {threads} &>{log}
         """

# index a pacbio bam file
rule index:
    output: "{basename}.bam.pbi"
    input: "{basename}.bam"
    conda: "conda/clustselect.yaml"
    envmodules:
        "bioinfo-tools",
        "SMRT/7.0.1"
    shell: "pbindex {input}"

# for each of the swarm clusters, make a BAM file containing the source CCS
checkpoint clustselect:
    output: directory("process/{clustype}/{seqrun}_{type}")
    input:
        clust="process/{seqrun}.ccs.{clustype}",
        uc=   "process/{seqrun}.ccs.derep.uc",
        bam=  expand("process/{movie}.{{type}}.bam", movie = moviefiles),
        pbi= expand("process/{movie}.{{type}}.bam.pbi", movie = moviefiles),
        script="scripts/clustselect.sh"
    log: "logs/{clustype}select_{seqrun}_{type}.log"
    threads: maxthreads
    conda: "conda/clustselect.yaml"
    #bamsieve in the env module gives mysterious errors
    #envmodules:
    #    "bioinfo-tools",
    #    "samtools",
    #    "SMRT/7.0.1",
    #    "gnuparallel/20180822"
    shell:
        """
	[ -d {output} ] || mkdir -p {output}
        cat {input.clust} |
          {{ parallel --pipe -N1 -j {threads} {input.script} {{#}} {input.uc} process .{wildcards.type} {output} {wildcards.clustype}; }} &>{log}
        """

# get consensus of CCS reads from each cluster
rule cluster_consensus:
    output: "process/{seqrun}.{clustype}.cons.fasta"
    input:
        clust="process/{seqrun}.ccs.{clustype}",
        uc=   "process/{seqrun}.ccs.derep.uc",
        fastq=  "process/{seqrun}.ccs.trimmed.fastq.gz",
        script="scripts/clust_consensus.sh",
        c3s= "bin/c3s"
    log: "logs/{clustype}_consensus_{seqrun}.log"
    threads: maxthreads
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1",
        "gnuparallel/20180822"
    shell:
        """
        cat {input.clust} |
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
    output: "process/{seqrun}_{type}.{clustype}.table"
    input:
        bamdir = "process/{clustype}/{seqrun}_{type}",
        script = "scripts/clust_table.sh",
        samples = "process/sample.tsv"
    log: "logs/{clustype}_table_{seqrun}_{type}.log"
    threads: maxthreads
    conda: "conda/samtools_parallel.yaml"
    envmodules:
        "bioinfo-tools",
        "samtools/1.10",
        "gnuparallel/20180822"
    shell:
        """
        ls -1 {input.bamdir}/{wildcards.clustype}_*.bam | {{ parallel -j {threads} {input.script} {{}} {input.samples}; }} >{output} 2>{log}
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

def reportfiles(wildcards):
    swarmdir = checkpoints.clustselect.get(**wildcards).output[0]
    clusters = glob_wildcards(os.path.join(swarmdir, "swarm_{cluster}.bam"))
    return expand("{swarmdir}/swarm_{cluster}.report.csv", swarmdir = swarmdir, cluster = clusters.cluster)

def laafastq(wildcards):
    swarmdir = checkpoints.clustselect.get(**wildcards).output[0]
    clusters = glob_wildcards(os.path.join(swarmdir, "swarm_{cluster}.bam"))
    return expand("{swarmdir}/swarm_{cluster}.laa.fastq.gz", swarmdir = swarmdir, cluster = clusters.cluster)

def swarmfiles(wildcards):
    swarmdir = checkpoints.clustselect.get(**wildcards).output[0]
    clusters = glob_wildcards(os.path.join(swarmdir, "swarm_{cluster}.bam"))
    return expand("{swarmdir}/swarm_{cluster}.{ext}", swarmdir = swarmdir, cluster = clusters.cluster,
                  ext = ["laa.fastq.gz", "junk.fastq.gz", "report.csv", "subreads.csv", "pcr.csv"])

rule laa_table:
    output: "process/{seqrun}_{type}.swarm.otutab"
    input:
        expand("process/{movie}.subreads.demux.sort.lima.report", movie= moviefiles),
        reportfiles,
        samplekey = "process/sample.tsv",
        script = "scripts/lima_map.sh"
    envmodules:
        "bioinfo-tools",
        "gnuparallel/20180822",
        "gawk/4.1.4"
    conda: "conda/gawk_parallel.yaml"
    threads: maxthreads
    log: "logs/laa_table_{seqrun}_{type}.log"
    shell:
        """
        ls -1 process/swarm/{wildcards.seqrun}_{wildcards.type}/*.subreads.csv |
        {{ parallel -j {threads} {input.script} "process" ".subreads.demux.sort" {input.samplekey} {{}}; }} >{output} 2>{log}
        """

rule laa_fastq:
    output: "process/{seqrun}_{type}.swarm.laa.fastq.gz"
    input: laafastq
    shell:
        """
        for f in process/swarm/{wildcards.seqrun}_{wildcards.type}/swarm_*.laa.fastq.gz; do
            zcat $f |
            sed 's/^@Barcode/@'"$(basename $f .laa.fastq.gz)"' Barcode/'
        done |
        gzip -c - >{output}
        """

rule laa_select:
    output: "process/{seqrun}_{type}.swarm.laa.select.fastq.gz"
    input:
        fastq = "process/{seqrun}_{type}.swarm.laa.fastq.gz",
        otutab = "process/{seqrun}_{type}.swarm.otutab"
    conda: "conda/vsearch.yaml"
    envmodules:
        "bioinfo-tools",
        "vsearch/2.14.1"
    threads: maxthreads
    log: "logs/laa_select_{seqrun}_{type}.log"
    shell:
        """
        temp=$(mktemp)
        trap 'rm $temp' EXIT
        awk '!x[$2 " " $3]++ {{print $2 " " $3}}' {input.otutab} > $temp
        vsearch --fastx_getseqs {input.fastq} --labels $temp --notrunclabels --threads {threads} --fastqout - 2>{log} |
        gzip -c - >{output}
        """

rule all_laa:
    output: touch("process/all_laa_{seqrun}_{type}")
    input: swarmfiles

rule r_targets:
    output: touch(".targets")
    input:
        expand(["process/pb_363.{cluster_type}.cons.fasta", "process/pb_363_ccs.{cluster_type}.table"],
            cluster_type = ["swarm", "vclust"]),
        "processReads/ampliseq/feature-table.tsv",
        expand("reference/{reference}.sintax.fasta.gz",
            reference = ["unite_single.ITS", "rdp_train.LSU", "silva_nr99.LSU"]),
        "reference/constraints.fasta",
        "processReads/ampliseq/qiime2_ASV_table.tsv",
        "start_files/meta.txt"
    conda: "conda/kungsangen.yaml"
    threads: maxthreads
    log: "logs/targets.log"
    shell:
        """
        R -e 'library(targets); tar_make(); stopifnot(all(is.na(tar_meta()$error)))' 2>1 | tee {log}
        """

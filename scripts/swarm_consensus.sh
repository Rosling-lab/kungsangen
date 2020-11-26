#!/usr/bin/env bash
# seqrun n uc ccs
seqrun="$1"
shift
n="$1"
shift
uc="$1"
shift
ccs="$1"
shift
conslog="logs/consswarm_$seqrun_$n.log"
echo "conslog=$conslog" >"$conslog"
echo "seqrun=$seqrun" >>"$conslog"
echo "n=$n" >>"$conslog"
echo "uc=$uc" >>"$conslog"
echo "ccs=$ccs" >>"$conslog"
echo "mapping swarm to reads" >>"$conslog"
templabel=$(mktemp)
tempfastq=$(mktemp --suffix ".fastq")
tempfasta=$(mktemp --suffix ".fasta")
trap 'rm $templabel $tempfastq $tempfasta' EXIT
tr " " "\\n" 2>>"$conslog" |
# remove abundances
sed 's/ccs;size=[0-9]*$/ccs/' 2>>"$conslog" |
# find all lines in the clustering file that contain one of these sequence ids
grep -F -f - "$uc" 2>>"$conslog" |
# S and H represent the first and subsequent sequences in the cluster
# column 9 is the query sequence
awk '/^[SH]/{sub(/\/ccs$/, "", $9); print $9}' >templabel 2>>"$conslog" &&
echo "extracting $(wc -l $templabel) reads" >>"$conslog" &&
# extract the sequences from the ccs file
vsearch --fastx_getseqs "$ccs" \
  --labels templabel\
  --fastqout "$tempfastq" &>>"$conslog" &&
echo "finding consensus" >>"$conslog" &&
bin/c3s "$tempfastq" "$tempfasta" &>>"$conslog" &&
sed 's/>consensus/&'$n'/' "$tempfasta" 2>>"$conslog"

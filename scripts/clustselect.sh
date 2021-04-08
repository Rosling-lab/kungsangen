#!/usr/bin/env bash
# helper script to extract given sequences from bam files
# the names of the sequences to select are taken from stdin
# command line parameters are:
# n uc bamdir bamsuffix outdir bamprefix ...
# n - integer index of this cluster
# uc - path to the .uc (uclust format) file from dereplicating the source sequences prior to clustering
# bamdir - directory where the input .bam files are found
# bamsuffix - a suffix to add to movie names (before .bam) to get the input .bam file names
# outdir - directory to put the output cluster files
# bamprefix - prefix to add to the output cluster files
# in summary, the input files are ${bamdir}/${moviename}${bamsuffix}.bam,
# where "moviename" is given by the sequence names up to the first "/",
# and the output file will be ${outdir}/${bamprefix}_${n}.bam (with $n formatted to 5 digits)
printf -v n "%05d" $1
echo "cluster $n"
shift
uc="$1"
shift
bamdir="$1"
shift
bamsuffix="$1"
shift
outdir="$1"
shift
bamprefix="$1"
shift
temp=$(mktemp -d swarmtempXXXXXXXX)
# put each sequence id on its own line
tr " " "\n" |
# remove abundances
sed -n 's/ccs;size=[0-9]*$/ccs/p' |
# find all lines in the clustering file that contain one of these sequence ids
grep -F -f - "$uc" |
# S and H represent the first and subsequent sequences in the cluster
# column 9 is the query sequence
awk '/^[SH]/{sub(/\/ccs$/, "", $9); print $9}' |
# split into a separate file for each movie
awk -F/ '{out= "'"$temp"'/" $1 "_'"$n"'.zmw"; print $2 > out}'

# extract matching ZMWs from each movie
allbams=""
for f in $(ls $temp/m*_$n.zmw); do
  b=$(basename -s "_$n.zmw" $f)
  newbam="$temp/$b.bam"
  bamsieve --whitelist $f "$bamdir/$b$bamsuffix.bam" "$newbam"
  allbams="$allbams $newbam"
done

# put everything together in one bam
[ -d "$outdir" ] || mkdir -p "$outdir"
pbmerge -o "$outdir/${bamprefix}_$n.bam" $allbams
rm -r "$temp"

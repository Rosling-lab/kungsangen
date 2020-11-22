#!/usr/bin/env bash
# n uc bamdir bamsuffix outdir ...
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
temp=$(mktemp -d swarmtempXXXXXXXX)
# put each sequence id on its own line
tr " " "\n" |
# remove abundances
sed 's/ccs;size=[0-9]*$/ccs/' |
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
pbmerge -o "$outdir/swarm_$n.bam" $allbams
rm -r "$temp"

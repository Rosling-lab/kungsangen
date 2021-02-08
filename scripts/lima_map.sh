#!/usr/bin/env bash
# keydir keyext barcodes files...
keydir=$1
shift
keyext=$1
shift
barcodes=$1
shift
temp=$(mktemp)
alltemp=$temp
trap 'rm $alltemp' EXIT
gawk -F "[/,]" -M -v "PREC=oct" -v "SUBSEP=/" '
  BEGINFILE {
    n=split(FILENAME, b, "[/.]")
    f=b[n-2]
  }
  FNR == 1 {
    for (i=2; i<=NF; i++)
      seq[i-1]=$i
    nclust=NF-1
    next
  }
  (! ($1, $2) in zmws) {
    zmws[$1,$2]=TRUE
    for (i = 1; i<=nclust; i++)
      p[$1,$2][i]=1
  }
  {
    for (i = 4; i <= NF; i++) {
      p[$1,$2][i-3] *= $i == 0 ? 1e-7 : $i
    }
  }
  END {
    for (z in p) {
      max=0
      total=0
      for (i in p[z]) {
        total+=p[z][i]
        if (p[z][i]>max) {
          max=p[z][i]
          maxi=i
        }
      }
      post=max/total
      if (max/total > 0.9)
        print z, f, seq[maxi], max/total
    }
  }' "$@" >$temp

movies=$(gawk -F / '{print $1}' $temp | sort -u) 

for m in $movies;
  do
  grep "^$m" $temp |
  sort -k 1b,1 |
  join - ${keydir}/${m}${keyext}.lima.report -o '1.2 1.3 2.41 2.42'
done |
  gawk '{print $3 "," $4, $1, $2}' |
  sort |
  uniq -c |
  gawk '{print $2, $3, $4, $1}' |
  join - $barcodes -o "2.2 1.2 1.3 1.4"
  

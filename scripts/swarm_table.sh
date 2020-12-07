#!/usr/bin/env bash
swarm=$(basename $1 .bam)
samtools view $1 |
  sed -n -r 's%(m[^/]+/[0-9]+).*bc:B:S,(1*[0-9],1[0-9]).*%\1 \2%p' |
  uniq |
  awk '{print $2}' |
  sort |
  uniq -c |
  awk '{print "'"${swarm}"'", $2, $1}'

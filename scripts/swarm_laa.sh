#!/usr/bin/env bash
# uc bam ...
uc="$1"
shift
bam="$1"
shift
temp=$(mktemp swarmtempXXXXXXXX)
tr " " "\n" |
sed 's/ccs;size=[0-9]*$//' |
grep -F -f - "$uc"

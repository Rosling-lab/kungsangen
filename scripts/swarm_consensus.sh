#!/usr/bin/env bash
# fasta swarmdir n logprefix
conslog="logs/consswarm_$4_$3.log"
echo "conslog=$conslog" >"$conslog"
swarmfasta="$2/swarm_$3.fasta"
echo "swarmfasta=$swarmfasta" >>"$conslog"
extractlog="logs/extractswarm_$4_$3.log"
echo "extractlog=$extractlog" >>"$conslog"
samplelog="logs/sampleswarm_$4_$3.log"
echo "samplelog=$samplelog" >>"$conslog"
alignlog="logs/alignswarm_$4_$3.log"
echo "alignlog=$alignlog" >>"$conslog"
alignfasta="$2/swarm_$3.aln.fasta"
echo "alignfasta=$alignfasta" >>"$conslog"
samplename="swarm$3"
echo "samplename=$samplename" >>"$conslog"
consfasta="$2/swarm_$3.cons.fasta"
echo "consfasta=$consfasta" >>"$conslog"
echo "extracting input swarm" >>"$conslog"
tr " " "\\n" |
vsearch --fastx_getseqs "$1" \
  --labels -\
  --fastaout "$swarmfasta" &>"$extractlog" &&
echo "counting sequences" >>"$conslog" &&
n=$(awk '/^>/{sub(/.+;size=/, "", $1); s+=$1}; END{print s}' "$swarmfasta") &&
echo "$n total sequences" >>"$conslog" &&
[ $n -eq 1 ] &&
sed '/^>/c >'"$samplename" "$swarmfasta" >"$consfasta" ||
{
  {
    [ $n -ge 100 ] &&
    vsearch --fastx_subsample "$swarmfasta"\
      --sizein\
      --sample_size 100\
      --fastaout - 2>"$samplelog" ||
    vsearch --rereplicate "$swarmfasta"\
      --output - 2>"$samplelog"
  } |
  muscle -maxiters 1 -diags -gapopen -0.5 -out "$alignfasta" 2>"$alignlog" &&
  cat "$alignfasta" |
  cons -filter -sformat fasta -osformat fasta -name "$samplename" 2>>"$conslog" |
  sed -n '/^>/{p;n;h;d}; H; ${x;s/[n\n]//g;p}' >"$consfasta" 2>>"$conslog"
}
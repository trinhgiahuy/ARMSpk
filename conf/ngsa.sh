#!/bin/bash

export APPDIR="./NGSAnalyzer"
export BINARY="./bin/workflow"
export INPUTDIR="./ngsa_mini_input"
export INPUT="$INPUTDIR/bwa_db/reference.fa $INPUTDIR/seq_contig.md $INPUTDIR/reference.fa $INPUTDIR/reference.fa.fai $INPUTDIR/00-read-rank"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6 2|48 2|24 2|12 2|6 4|12 4|6 6|8 6|4 6|2 12|2 12|4 12|3 12|4 12|5 12|6"
	export BESTCONF=""
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

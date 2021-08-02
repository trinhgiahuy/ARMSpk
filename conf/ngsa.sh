#!/bin/bash

export APPDIR="./NGSAnalyzer"
export BINARY="./bin/workflow"
export INPUTDIR="./ngsa_mini_input"
export INPUT="$INPUTDIR/bwa_db/reference.fa $INPUTDIR/seq_contig.md $INPUTDIR/reference.fa $INPUTDIR/reference.fa.fai $INPUTDIR/00-read-rank"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="140m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export MAXTIME="100m"
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 2|6 2|12 2|24
			 4|1 4|2 4|4 4|6 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1 48|2
			 96|1"
	export BESTCONF="12|4"
	export SCALCONF="12|84" #only up to 12 MPI ranks supported for the input... 32|32 128|8"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64 1|128 1|192 1|256
			 4|16 4|32 4|48 4|64
			 16|4 16|8 16|12 16|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6"
	export BESTCONF="4|32"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|64 1|72 1|128 1|144 1|192 1|256 1|288
			 4|18 4|36 4|54 4|72
			 16|4 16|8 16|12 16|16
			 18|4 18|8 18|12 18|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 72|1 72|2 72|4"
	export BESTCONF="4|18"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|36 1|48
			 2|6 2|8 2|12 2|16 2|24
			 4|4 4|6 4|8 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1 48|2
			 96|1"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF=""
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF=""
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF=""
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF=""
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF=""
	export NumRunsGEM5=1
fi

#!/bin/bash

export APPDIR="./MODYLAS/data/wat222"
export BINARY="../../src/modylas_mini"
export INPUT="./wat222"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="8|1 8|2 8|3 8|4 8|6 8|8 8|10
			 16|1 16|2 16|3 16|4
			 32|1 32|2
			 64|1"
	export BESTCONF="8|6" #openmpi maps this config oddly but fallback is as fast as before with intel mpi"16|3"
	export SCALCONF="8|128 32|32 128|8"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="8|4 8|8 8|16 8|32
			 16|4 16|8 16|12 16|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 128|1 128|2 128|3
			 256|1 256|2"
	export BESTCONF="64|4"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="8|4 8|8 8|16 8|32
			 16|4 16|8 16|12 16|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 128|1 128|2 128|3
			 256|1 256|2"
	export BESTCONF="64|4"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="8|1 8|2 8|3 8|4 8|6
			 16|1 16|2 16|3 16|4
			 32|1 32|2"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="16|3"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="16|3"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="16|3"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="16|3"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF=""	#no single rank option
	export NumRunsGEM5=1
fi

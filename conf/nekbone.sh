#!/bin/bash

export APPDIR="./Nekbone/test/nek_mgrid"
export BINARY="./nekbone"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export ielN=8640	#2^6*3^3*5 to split nicely -> mcmodel=large needed above 10k
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="6|1 12|1 16|1 24|1 32|1 48|1 96|1"
	export BESTCONF="48|1" # no fucking clue and so sick of computers "96|1"
	export SCALCONF="48|21 128|8"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="16|1 32|1 64|1 96|1 128|1 192|1 256|1"
	export BESTCONF="128|1"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="64|1 72|1 96|1 128|1 144|1 192|1 256|1 288|1"
	export BESTCONF="72|1"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|36 1|48
			 2|4 2|8 2|12 2|16 2|24
			 4|4 4|6 4|8 4|12
			 6|2 6|4 6|8
			 12|1 12|2 12|4
			 16|1 16|2 16|4
			 24|1 24|2
			 32|1 32|2
			 48|1"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="48|1"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="48|1"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="48|1"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="48|1"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|1"
	export NumRunsGEM5=1
fi

#!/bin/bash

export APPDIR="./BabelStream"
export BINARYS="./omp-stream_2 ./omp-stream_8" #./omp-stream_14"
export INPUT="-s SIZE -n 10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="20m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export MAXTIME="100m"
	export TESTCONF="1|2 1|6 1|12 1|18 1|24"
	export BESTCONF="1|12"
	export SCALCONF="1|256 1|1024"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="1|8 1|16 1|32 1|64 1|96 1|128"
	export BESTCONF="1|64"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|9 1|18 1|36 1|72 1|108 1|144"
	export BESTCONF="1|144"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|32 1|36 1|48"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="1|36"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="1|36"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="1|24"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="1|32"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|36"
	export NumRunsGEM5=1
fi

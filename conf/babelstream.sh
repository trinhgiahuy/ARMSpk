#!/bin/bash

export APPDIR="./BabelStream"
export BINARYS="./omp-stream_2 ./omp-stream_14"
export INPUT="-s SIZE -n 10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="2m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export MAXTIME="10m"
	export TESTCONF="2 6 12 18 24"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="8 16 32 64 96 128"
	export BESTCONF="64"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="9 18 36 72 108 144"
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

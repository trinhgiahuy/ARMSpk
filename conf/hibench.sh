#!/bin/bash

export APPDIR="./HiBench"
export BINARYS="./bin/run_all.sh"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="60m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1"
	export BESTCONF="1"
	export SCALCONF=""
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1"
	export BESTCONF="1"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="1"
	export BESTCONF="1"
else
	echo "Unsupported host"
	exit
fi

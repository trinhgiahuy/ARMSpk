#!/bin/bash

export APPDIR="./Nekbone/test/nek_mgrid"
export BINARY="./nekbone"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export ielN=1024
export MAXTIME="1m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="no" #"yes"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="6|1 12|1 16|1 24|1 32|1 48|1 96|1"
	export BESTCONF="24|1"
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="16|1 32|1 64|1 96|1 128|1 192|1 256|1"
	export BESTCONF="64|1"
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

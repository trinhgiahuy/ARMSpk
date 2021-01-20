#!/bin/bash

export APPDIR="./Nekbone/test/nek_mgrid"
export BINARY="./nekbone"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export ielN=13824
export MAXTIME="1m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="6|1 12|1 16|1 24|1 32|1 48|1 96|1"
	export BESTCONF="48|1" # no fucking clue and so sick of computers "96|1"
	export SCALCONF="48|21 128|8"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="16|1 32|1 64|1 96|1 128|1 192|1 256|1"
	export BESTCONF="128|1"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="64|1 72|1 96|1 128|1 144|1 192|1 256|1 288|1"
	export BESTCONF="72|1"
else
	echo "Unsupported host"
	exit
fi

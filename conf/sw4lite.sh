#!/bin/bash

export APPDIR="./SW4lite"
export INPUT="./tests/pointsource/pointsource.in"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export BINARY="./optimize_mp_kiev/sw4lite"
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF=""
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export BINARY="./optimize_mp_mill/sw4lite"
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export BINARY="./optimize_mp_mill/sw4lite"
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

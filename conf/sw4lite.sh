#!/bin/bash

export APPDIR="./SW4lite"
export INPUT="./tests/pointsource/pointsource.in"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export BINARY="./optimize_mp_kiev/sw4lite"
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 2|6 2|12 2|24
			 4|1 4|2 4|4 4|6 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1
			 96|1"
	export BESTCONF="24|2"
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

#!/bin/bash

export APPDIR="./MACSio"
export BINARY="./macsio/macsio"
export MAXNDPP=$((2*2*2*2*2*3*3))
export INPUT="--units_prefix_system decimal --num_dumps NDPP"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="2m"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|1 4|1 6|1 12|1 24|1 32|1 48|1 96|1"
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

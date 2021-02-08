#!/bin/bash

export APPDIR="./MACSio"
export BINARY="./macsio/macsio"
export MAXNDPP=1024  #increase by ~3x to support 1k MPI threads; $((2*2*2*2*2*3*3))
# but decrease size of dump by ~3x from 80k to 20k
export INPUT="--part_size 20000 --units_prefix_system decimal --num_dumps NDPP"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="2m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|1 4|1 6|1 12|1 24|1 32|1 48|1 96|1"
	export BESTCONF="4|1"
	export SCALCONF="256|1 1024|1"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="16|1 32|1 64|1 96|1 128|1 192|1 256|1"
	export BESTCONF="64|1"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="16|1 18|1 32|1 64|1 72|1 96|1 128|1 144|1 192|1 256|1 288|1"
	export BESTCONF="64|1"
else
	echo "Unsupported host"
	exit
fi

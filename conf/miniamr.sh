#!/bin/bash

export APPDIR="./MiniAMR/ref"
export BINARY="./ma.x"
export INPUT="--num_refine 4 --max_blocks 9000 --npx PX --npy PY --npz PZ --nx 2 --ny 2 --nz 2 --num_objects 1 --object 2 0 -1.71 -1.71 -1.71 0.04 0.04 0.04 1.7 1.7 1.7 0.0 0.0 0.0 --num_tsteps 100 --checksum_freq 1"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="7m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="2|1|2|1|1 4|1|2|2|1 6|1|3|2|1 12|1|3|2|2 24|1|4|3|2 32|1|4|4|2 48|1|4|4|3 96|1|6|4|4"
	export BESTCONF="96|1|6|4|4"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="32|1|4|4|2 48|1|4|4|3 64|1|4|4|4 80|1|5|4|4 96|1|6|4|4 128|1|8|4|4 192|1|8|6|4 256|1|8|8|4"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="32|1|4|4|2 36|1|6|3|2 48|1|4|4|3 64|1|4|4|4 72|1|6|4|3 80|1|5|4|4 96|1|6|4|4 128|1|8|4|4 144|1|6|6|4 192|1|8|6|4 256|1|8|8|4 288|1|8|6|6"
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

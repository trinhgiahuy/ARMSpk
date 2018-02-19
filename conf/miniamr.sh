#!/bin/bash

export APPDIR="./MiniAMR/ref"
export BINARY="./ma.x"
export INPUT="--num_refine 4 --max_blocks 9000 --npx PX --npy PY --npz PZ --nx 8 --ny 8 --nz 8 --num_objects 1 --object 2 0 -1.71 -1.71 -1.71 0.04 0.04 0.04 1.7 1.7 1.7 0.0 0.0 0.0 --num_tsteps 100 --checksum_freq 1"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="6|1|3|2|1 12|1|3|2|2 24|1|4|3|2 32|1|4|4|2 48|1|4|4|3 96|1|6|4|4"
	export BESTCONF="48|1|4|4|3"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

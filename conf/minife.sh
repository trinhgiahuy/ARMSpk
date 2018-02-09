#!/bin/bash

export APPDIR="./MiniFE"
export BINARYS="./mkl/src/miniFE.x ./openmp-opt-knl/src/miniFE.x ./openmp-opt/src/miniFE.x"
export INPUT="-nx 128 -ny 128"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1 96|1"
	BBINARY="./openmp-opt/src/miniFE.x"
	export BESTCONF="24|1"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

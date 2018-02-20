#!/bin/bash

export APPDIR="./MiniFE"
export BINARYS="./mkl/src/miniFE.x ./openmp-opt-knl/src/miniFE.x ./openmp-opt/src/miniFE.x"
export BBINARY="./openmp-opt/src/miniFE.x"
export INPUT="-nx 128 -ny 128"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1 96|1"
	export BESTCONF=""
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

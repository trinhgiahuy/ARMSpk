#!/bin/bash

export APPDIR="./MACSio"
export BINARY="./macsio/macsio"
export MAXNDPP=1024  #increase by ~3x to support 1k MPI threads; $((2*2*2*2*2*3*3))
# but decrease size of dump by ~3x from 80k to 20k
export INPUT="--part_size 20000 --units_prefix_system decimal --num_dumps NDPP"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="20m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|1 4|1 6|1 12|1 24|1 32|1 48|1 96|1"
	export BESTCONF="4|1"
	export SCALCONF="256|1 1024|1"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="16|1 32|1 64|1 96|1 128|1 192|1 256|1"
	export BESTCONF="64|1"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="16|1 18|1 32|1 64|1 72|1 96|1 128|1 144|1 192|1 256|1 288|1"
	export BESTCONF="64|1"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|1 4|1 6|1 12|1 24|1 36|1 48|1"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="48|1"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="48|1"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="48|1"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="48|1"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|1"
	export NumRunsGEM5=1
fi

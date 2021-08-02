#!/bin/bash

export APPDIR="./DLproxy/benchmarks/conv_gemm"
export BINARY="./main"
export INPUT="FP32 32 3 224 3 32 1 10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
        # on "normal" Xeon
        export TESTCONF="1|4 1|8 1|12 1|24 1|32 1|36 1|48"
        export BESTCONF="1|24"
elif [ -n "${IKNLHOST}" ]; then
        # on one of the Phi (knl)
        export TESTCONF=""
        export BESTCONF=""
elif [ -n "${IKNMHOST}" ]; then
        # on one of the Phi (knm)
        export TESTCONF=""
        export BESTCONF=""
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|24 1|32 1|36 1|48"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="1|48"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="1|48"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="1|48"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="1|48"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|48"
	export NumRunsGEM5=1
fi

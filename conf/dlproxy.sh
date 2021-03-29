#!/bin/bash

export APPDIR="./DLproxy/benchmarks/conv_gemm"
export BINARY="./main"
export INPUT="FP32 32 3 224 3 32 1 10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="1m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
        # on "normal" Xeon
        export TESTCONF="24"
        export BESTCONF="24"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
        # on one of the Phi (knl)
        export TESTCONF=""
        export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
        # on one of the Phi (knm)
        export TESTCONF=""
        export BESTCONF=""
else
        echo "Unsupported host"
        exit
fi

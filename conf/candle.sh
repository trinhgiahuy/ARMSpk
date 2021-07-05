#!/bin/bash

export APPDIR="./CANDLE"
export BINARYS="p1b1_baseline_keras2.py" # p1b2_baseline_keras2.py p1b3_baseline_keras2.py p2b1_baseline_keras2.py p2b2_baseline_keras2.py p3b1_baseline_keras2.py p3b2_baseline_keras2.py"
export INPUT=""
export PATH=$ROOTDIR/dep/anaconda2/bin:$PATH
export MKL_THREADING_LAYER=GNU
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="no"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export MAXTIME="10m"
	export TESTCONF="1|2 1|6 1|12 1|24 1|48"
	export BESTCONF="1|12"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="1|32 1|64 1|128 1|192 1|256"
	export BESTCONF="1|32"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|32 1|36 1|64 1|72 1|144"
	export BESTCONF="1|64"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF=""
fi

#!/bin/bash

export APPDIR="./CANDLE"
export BINARYS="p1b1_baseline_keras2.py" # p1b2_baseline_keras2.py p1b3_baseline_keras2.py p2b1_baseline_keras2.py p2b2_baseline_keras2.py p3b1_baseline_keras2.py p3b2_baseline_keras2.py"
export INPUT=""
export PATH=$ROOTDIR/dep/anaconda2/bin:$PATH
export MKL_THREADING_LAYER=GNU
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="40m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1 2 6 12 24 48"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1 32 64 128 192 256"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="1 36 72 144 216 288"
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

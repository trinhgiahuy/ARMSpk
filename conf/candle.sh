#!/bin/bash

export APPDIR="./CANDLE"
export BINARYS="p1b1_baseline_keras2.py" # p1b2_baseline_keras2.py p1b3_baseline_keras2.py p2b1_baseline_keras2.py p2b2_baseline_keras2.py p3b1_baseline_keras2.py p3b2_baseline_keras2.py"
export PATH=$ROOTDIR/dep/anaconda2/bin:$PATH
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1"
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

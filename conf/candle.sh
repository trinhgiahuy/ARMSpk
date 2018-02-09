#!/bin/bash

export APPDIR="./CANDLE"
export BINARYS="Pilot1/P1B1/p1b1_baseline_keras2.py" #Pilot1/P1B2/p1b2_baseline_keras2.py Pilot1/P1B3/p1b3_baseline_keras2.py Pilot2/P2B1/p2b1_baseline_keras2.py Pilot2/P2B2/p2b2_baseline_keras2.py Pilot3/P3B1/p3b1_baseline_keras2.py Pilot3/P3B2/p3b2_baseline_keras2.py"
export PATH=$HOME/anaconda2/bin:$PATH

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1"
	export BESTCONF="1"
else
	# on one of the Phi
	export TESTCONF="1"
	export BESTCONF="1"
fi

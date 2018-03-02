#!/bin/bash

export APPDIR="./HPCG/build"
export BINARY="./bin/xhpcg"
export MAXXYZ=$((2*2*2*3*3*5))
export INPUT="--nx=NX --ny=NY --nz=NZ"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|12 1|24 1|32 1|48
			 2|6 2|12 2|18 2|24
			 12|1 12|2 12|3 12|4
			 24|1 24|2"
	export BESTCONF="2|24"
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

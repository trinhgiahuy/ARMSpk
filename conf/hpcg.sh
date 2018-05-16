#!/bin/bash

export APPDIR="./HPCG/build"
export BINARY="./bin/xhpcg"
export MAXXYZ=$((2*2*2*3*3*5))
export INPUT="--nx=NX --ny=NY --nz=NZ"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="no" #"yes"
export RUNPCM="yes"
export RUNVTUNE="no" #"yes"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|12 1|24 1|32 1|48
			 2|6 2|12 2|18 2|24
			 12|1 12|2 12|3 12|4
			 24|1 24|2"
	export BESTCONF="2|24"
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64 1|128 1|192 1|256
			 4|16 4|32 4|48 4|64
			 16|4 16|8 16|12 16|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 96|1 96|2 96|3
			 128|1 128|2
			 192|1
			 256|1"
	export BESTCONF="64|1"
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF="64|1" # temporary for testing freq scaling
else
	echo "Unsupported host"
	exit
fi

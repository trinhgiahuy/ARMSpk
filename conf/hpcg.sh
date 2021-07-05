#!/bin/bash

export APPDIR="./HPCG/build"
export BINARY="./bin/xhpcg"
#export MAXXYZ=$((2*2*2*3*3*5)) #doesnt work in flat mode
export MAXXYZ=$((2*2*2*3*5))
export INPUT="--nx=NX --ny=NY --nz=NZ"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|12 1|24 1|32 1|48
			 2|6 2|12 2|18 2|24
			 12|1 12|2 12|3 12|4
			 24|1 24|2"
	export BESTCONF="2|24"
	export SCALCONF="2|512 32|32 128|8"
elif [ -n "${IKNLHOST}" ]; then
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
	export BESTCONF="96|1"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|64 1|72 1|128 1|144 1|192 1|256 1|288
			 4|18 4|36 4|54 4|72
			 16|4 16|8 16|12 16|16
			 18|4 18|8 18|12 18|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 72|1 72|2 72|4
			 96|1 96|2 96|3
			 128|1 128|2
			 144|1 144|2
			 192|1
			 256|1
			 288|1"
	export BESTCONF="72|1"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|36 1|48
			 2|6 2|8 2|12 2|16 2|24
			 4|4 4|6 4|8 4|12
			 12|1 12|2 12|3 12|4
			 24|1 24|2
			 48|1"
	export BESTCONF=""
	export SCALCONF=""
fi

#!/bin/bash

export APPDIR="./Laghos"
export BINARY="./laghos"
#export INPUT="-p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8 -no-vis -pa"
export INPUT="-p 1 -pt 311 -m data/cube01_hex.mesh -rs 2 -tf 0.1 -no-vis -pa" #XXX: this only works for 24 mpi procs, but the non-cartesian code is broken AF with openmpi; also scale down runtime with tf 06 -> tf 0.1
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="30m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 2|6 2|12 2|24
			 4|1 4|2 4|4 4|6 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1
			 96|1"
	export BESTCONF="24|1"
	export SCALCONF="24|42 81|12" #partitioning is complex and with 311 only 3, 24, 81, ... (3*(x^3), x=1,2,3,...) are valid ###32|32 128|8"
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
	export BESTCONF="64|4"
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
	export BESTCONF="64|1"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|36 1|48
			 2|6 2|8 2|12 2|16 2|24
			 4|1 4|2 4|4 4|6 4|8 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="48|1"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="48|1"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="48|1"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="48|1"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|1"
	export NumRunsGEM5=1
fi

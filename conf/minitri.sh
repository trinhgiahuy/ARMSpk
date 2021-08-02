#!/bin/bash

export APPDIR="./MiniTri"
export BINARYMPI="./miniTri/linearAlgebra/MPI/miniTri.exe"
export BINARYOMP="./miniTri/linearAlgebra/openmp/miniTri.exe"
export INPUTMPI="./bcsstk30.mtx MM"
export INPUTOMP="./bcsstk30.mtx 16 OMPNT MM"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 4|1
			 6|1
			 12|1
			 24|1
			 32|1
			 48|1
			 96|1"
	export BESTCONF="1|48"
	export SCALCONF="1|1024 1024|1"
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64 1|128 1|192 1|256
			 16|1
			 32|1
			 64|1
			 96|1
			 128|1
			 192|1
			 256|1"
	export BESTCONF="1|128"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|64 1|72 1|128 1|144 1|192 1|256 1|288
			 64|1
			 72|1
			 96|1
			 128|1
			 144|1
			 192|1
			 256|1
			 288|1"
	export BESTCONF="1|128"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export MAXTIME="110m"
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|36 1|48
			 4|1
			 6|1
			 8|1
			 12|1
			 24|1
			 32|1
			 48|1"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="32|1"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="32|1"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="1|48"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="1|48"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|1"
	export NumRunsGEM5=1
fi

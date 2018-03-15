#!/bin/bash

export APPDIR="./MiniTri"
export BINARYMPI="./miniTri/linearAlgebra/MPI/miniTri.exe"
export BINARYOMP="./miniTri/linearAlgebra/openmp/miniTri.exe"
export INPUTMPI="./bcsstk30.mtx MM"
export INPUTOMP="./bcsstk30.mtx 16 OMPNT MM"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="1m"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"kiev"* ]]; then
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
elif [[ $HOSTNAME = *"lyon"* ]]; then
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
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

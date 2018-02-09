#!/bin/bash

export APPDIR="./MiniTri"
export BINARYMPI="./miniTri/linearAlgebra/MPI/miniTri.exe"
export BINARYOMP="./miniTri/linearAlgebra/openmp/miniTri.exe"
export INPUTMPI="./bcsstk30.mtx MM"
export INPUTOMP="./bcsstk30.mtx 16 OMPNT MM"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6 6|1 12|1 24|1 32|1 48|1 96|1"
	export BESTCONF="1|48"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

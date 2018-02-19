#!/bin/bash

export APPDIR="$1/NTChem"
export BINARY="./bin/rimp2.exe"
export INPUT=""
export NTCHEM_DIR=$APPDIR
export MODEL="h2o"
export DATA_DIR=${NTCHEM_DIR}/tests/${MODEL}

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 2|24 2|12 4|12 4|6 8|3 8|6 12|1 12|2 12|4 24|1 24|2 32|1 32|2 48|1 48|2"
	export BESTCONF=""
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

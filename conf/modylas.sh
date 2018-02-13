#!/bin/bash

export APPDIR="./MODYLAS/data/wat222"
export BINARY="../../src/modylas_mini"
export INPUT="./wat222"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="8|1 8|2 8|3 8|6 12|1 12|2 12|4 24|1 24|1 32|1 48|1"
	export BESTCONF="32|1"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

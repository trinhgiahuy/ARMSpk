#!/bin/bash

export APPDIR="./MODYLAS/data/wat222"
export BINARY="../../src/modylas_mini"
export INPUT="./wat222"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="8|1 8|2 8|3 8|4 8|6 16|1 16|2 16|3 16|4 32|1 32|2 64|1"
	export BESTCONF=""
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

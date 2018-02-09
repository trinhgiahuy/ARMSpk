#!/bin/bash

export APPDIR="./XSBench/src"
export BINARY="./XSBench"
export INPUT="-t OMPNT -s large -l 15000000 -G unionized"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF="2|24"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

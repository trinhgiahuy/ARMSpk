#!/bin/bash

export APPDIR="./AMG"
export BINARY="./test/amg"
export INPUT=""

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF="1|12"
else
	# on one of the Phi
	export TESTCONF="" #1|256 1|128 1|64 2|128 2|64 2|32 4|64 4|32 4|16"
	export BESTCONF=""
fi

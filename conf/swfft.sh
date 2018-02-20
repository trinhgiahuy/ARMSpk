#!/bin/bash

export APPDIR="./SWFFT"
export INPUT="32 128"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export BINARY="./build.xeon/TestFDfft"
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF=""
else
	# on one of the Phi
	export BINARY="./build.xmic/TestFDfft"
	export TESTCONF=""
	export BESTCONF=""
fi

#!/bin/bash

export APPDIR="./SW4lite"
export INPUT="./tests/pointsource/pointsource.in"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export BINARY="./optimize_mp_kiev/sw4lite"
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF="24|1"
else
	# on one of the Phi
	export BINARY="./optimize_mp_mill/sw4lite"
	export TESTCONF=""
	export BESTCONF=""
fi

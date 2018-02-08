#!/bin/bash

export APPDIR="./MACSio"
export BINARY="./macsio/macsio"
export INPUT=""

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|1 12|1 24|1 32|1 48|1"
	export BESTCONF="24|1"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

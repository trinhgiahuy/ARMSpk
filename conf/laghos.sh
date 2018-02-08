#!/bin/bash

export APPDIR="./Laghos"
export BINARY="./laghos"
export INPUT="-p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8 -no-vis -pa"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="6|1 8|1 12|2 16|1 24|1 32|1 48|1"
	export BESTCONF="24|1"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

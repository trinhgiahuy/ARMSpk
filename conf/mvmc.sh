#!/bin/bash

export APPDIR="./MVMC/job_middle"
export BINARY="./vmc.out"
export INPUT="./multiDir.def"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6 2|24 2|12 2|6 4|12 4|6 6|8 6|4 12|4 12|2 24|1 24|2 32|1 32|2 48|1 48|2"
	export BESTCONF="4|6"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi
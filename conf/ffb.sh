#!/bin/bash

export APPDIR="./FFB/bin/"
export BINARY="./ffb_mini"
export INPUT="DX DY DZ 48"

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="6|1|3|2|1 12|1|3|2|2 24|1|4|3|2 32|1|4|4|2 48|1|4|4|3 96|1|6|4|4 6|2|3|2|1 12|2|3|2|2 24|2|4|3|2 32|2|4|4|2 48|2|4|4|3 6|4|3|2|1 12|4|3|2|2 6|8|3|2|1 2|6|2|1|1 2|12|2|1|1 2|24|2|1|1 1|6|1|1|1 1|12|1|1|1 1|24|1|1|1 1|32|1|1|1 1|48|1|1|1 1|96|1|1|1"
	export BESTCONF="1|6|1|1|1"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

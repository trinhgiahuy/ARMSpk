#!/bin/bash

export APPDIR="./AMG"
export BINARY="./test/amg"
export MAXXYZ=$((2*2*2*2*2*3*3))
export INPUT="-problem 1 -P PX PY PZ -n NX NY NZ"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|6|1|1|1 1|12|1|1|1 1|24|1|1|1 1|32|1|1|1 1|48|1|1|1 1|96|1|1|1
	                 2|6|2|1|1 2|12|2|1|1 2|24|2|1|1
			 4|1|2|2|1 4|2|2|2|1 4|4|2|2|1 4|6|2|2|1 4|12|2|2|1
			 6|1|3|2|1 6|2|3|2|1 6|4|3|2|1 6|8|3|2|1
			 12|1|3|2|2 12|2|3|2|2 12|4|3|2|2
			 24|1|4|3|2 24|2|4|3|2
			 32|1|4|4|2 32|2|4|4|2
			 48|1|4|4|3
			 96|1|6|4|4"
	export BESTCONF=""
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

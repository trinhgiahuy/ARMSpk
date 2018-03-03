#!/bin/bash

export APPDIR="./FFB/bin/"
export BINARY="./ffb_mini"
export MAXDCZ=$((50*50*50))
export INPUT="PX PY PZ DCZ"
export NumRunsTEST=3
export NumRunsBEST=10

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
	export BESTCONF="24|1|4|3|2"
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64|1|1|1 1|128|1|1|1 1|192|1|1|1 1|256|1|1|1
			 4|16|2|2|1 4|32|2|2|1 4|48|2|2|1 4|64|2|2|1
			 16|4|4|2|2 16|8|4|2|2 16|12|4|2|2 16|16|4|2|2
			 32|2|4|4|2 32|4|4|4|2 32|6|4|4|2 32|8|4|4|2
			 64|1|4|4|4 64|2|4|4|4 64|4|6|4|4 64|6|6|4|4
			 96|1|6|4|4 96|2|6|4|4 96|3|6|4|4
			 128|1|8|4|4 128|2|8|4|4
			 192|1|8|6|4
			 256|1|8|8|4"
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

#!/bin/bash

export APPDIR="./FFVC"
export BINARY="./bin/ffvc_mini"
export INPUT="--scale=strong --size=144 --division=DXxDYxDZ"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="6|1|3|2|1 12|1|3|2|2 24|1|4|3|2 32|1|4|4|2 48|1|4|4|3 96|1|6|4|4 6|2|3|2|1 12|2|3|2|2 24|2|4|3|2 32|2|4|4|2 48|2|4|4|3 6|4|3|2|1 12|4|3|2|2 6|8|3|2|1 2|6|2|1|1 2|12|2|1|1 2|24|2|1|1 1|6|1|1|1 1|12|1|1|1 1|24|1|1|1 1|32|1|1|1 1|48|1|1|1 1|96|1|1|1"
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

#!/bin/bash

export APPDIR="./HPL/bin/Linux_Intel64"
export BINARY="./xhpl"
export HPLNS=$((336*192))
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|12|1|1 1|24|1|1 1|32|1|1 1|48|1|1
			 2|6|2|1 2|12|2|1 2|18|2|1 2|24|2|1
			 12|1|4|3 12|2|4|3 12|3|4|3 12|4|4|3
			 24|1|6|4 24|2|6|4"
	export BESTCONF=""
	export HPLNB="192"
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
	export HPLNB="336"
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
	export HPLNB="336"
else
	echo "Unsupported host"
	exit
fi

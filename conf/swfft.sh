#!/bin/bash

export APPDIR="./SWFFT"
export INPUT="32 128"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export BINARY="./build.xeon/TestFDfft"
	export TESTCONF="1|96 1|48 1|24 1|12 2|24 2|12 4|12 4|6 12|2 24|1 32|1 48|1"
	export BESTCONF=""
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export BINARY="./build.xmic/TestFDfft"
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export BINARY="./build.xmic/TestFDfft"
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

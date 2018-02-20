#!/bin/bash

export APPDIR="./MVMC"
export BINARY="../src/vmc.out"
export INPUT="./multiDir.def"
export PATH=$HOME/anaconda2/bin:$PATH

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6 2|24 2|12 2|6 4|12 4|6 6|8 6|4 12|4 12|2 24|1 24|2 32|1 32|2 48|1 48|2"
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

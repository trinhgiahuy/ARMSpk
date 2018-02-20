#!/bin/bash

export APPDIR="./NICAM/test/case/jablonowski"
export BINARY="./nhm_driver"
export NICAM_SYS=Linux64-intel-impi
export INPUT="./gl05rl00z40pe10"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="10|1 10|2 10|3 10|4 10|5 10|6"
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

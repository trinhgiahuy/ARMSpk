#!/bin/bash

export APPDIR="./NICAM/test/case/jablonowski"
export BINARY="./nhm_driver"
export NICAM_SYS=Linux64-intel-impi
export INPUT="./gl05rl00z40pe10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="2m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="10|1 10|2 10|3 10|4 10|5 10|6"
	export BESTCONF="10|6"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="10|5 10|6 10|7 10|10 10|15 10|20 10|25"
	export BESTCONF="10|15"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="10|5 10|6 10|7 10|10 10|14 10|16 10|20 10|28 10|30"
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

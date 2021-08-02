#!/bin/bash

export APPDIR="./NICAM/test/case/jablonowski"
export BINARY="./nhm_driver"
export NICAM_SYS=Linux64-intel-impi
export INPUT="./gl05rl00z40pe10"
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="20m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="10|1 10|2 10|3 10|4 10|5 10|6"
	export BESTCONF="10|1" # doesnt seem to use threads at all, except on K "10|6"
	export SCALCONF="" # doesnt support strong scaling AND no omp... facepalm
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="10|5 10|6 10|7 10|10 10|15 10|20 10|25"
	export BESTCONF="10|15"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="10|5 10|6 10|7 10|10 10|14 10|16 10|20 10|28 10|30"
	export BESTCONF="10|20"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="10|1 10|2 10|3 10|4 10|5 10|6"
	if   [[ "$1" = *"fujitrad"* ]];  then export BESTCONF="10|4"
	elif [[ "$1" = *"fujiclang"* ]]; then export BESTCONF="10|4"
	elif [[ "$1" = *"llvm12"* ]];    then export BESTCONF="10|4"
	elif [[ "$1" = *"gnu"* ]];       then export BESTCONF="10|4"
	fi
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF=""	#no single rank option
	export NumRunsGEM5=1
fi

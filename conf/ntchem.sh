#!/bin/bash

export APPDIR="$1/NTChem"
export BINARY="./bin/rimp2.exe"
export INPUT=""
export NTCHEM_DIR=$APPDIR
export MODEL="h2o"
export DATA_DIR=${NTCHEM_DIR}/tests/${MODEL}
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="1m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="no" #"yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 2|6 2|12 2|24
			 4|1 4|2 4|4 4|6 4|12
			 6|1 6|2 6|4 6|8
			 12|1 12|2 12|4
			 24|1 24|2
			 32|1 32|2
			 48|1 48|2
			 96|1"
	export BESTCONF="24|1"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64 1|128 1|192 1|256
			 4|16 4|32 4|48 4|64
			 16|4 16|8 16|12 16|16
			 32|2 32|4 32|6 32|8
			 64|1 64|2 64|4 64|6
			 96|1 96|2 96|3
			 128|1 128|2
			 192|1
			 256|1"
	export BESTCONF="16|8"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

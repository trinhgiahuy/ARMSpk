#!/bin/bash

export APPDIR="./HPL/bin/Linux_Intel64"
export BINARY="./xhpl"
export HPLNS=$((336*192))
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|12|1|1 1|24|1|1 1|32|1|1 1|48|1|1
			 2|6|2|1 2|12|2|1 2|18|2|1 2|24|2|1
			 12|1|4|3 12|2|4|3 12|3|4|3 12|4|4|3
			 24|1|6|4 24|2|6|4"
	export BESTCONF="24|1|6|4"
	export HPLNB="192"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|48|1|1 1|64|1|1 1|96|1|1 1|128|1|1 1|256|1|1
			 2|32|2|1 2|64|2|1 2|96|2|1 2|128|2|1
			 4|8|2|2 4|16|2|2 4|32|2|2 4|64|2|2
			 8|4|4|2 8|8|4|2 8|16|4|2
			 16|4|4|4 16|6|4|4 16|8|4|4
			 32|1|8|4 32|2|8|4 32|3|8|4
			 64|1|8|8 64|2|8|8
			 128|1|16|8
			 256|1|16|16"
	export BESTCONF=""
	export HPLNB="336"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF="1|64|1|1 1|72|1|1 1|128|1|1 1|144|1|1 1|192|1|1 1|256|1|1 1|288|1|1
			 4|18|2|2 4|36|2|2 4|54|2|2 4|72|2|2
			 16|4|4|4 16|8|4|4 16|12|4|4 16|16|4|4
			 18|4|6|3 18|8|6|3 18|12|6|3 18|16|6|3
			 32|2|8|4 32|4|8|4 32|6|8|4 32|8|8|4
			 64|1|8|8 64|2|8|8 64|4|8|8 64|6|8|8
			 72|1|9|8 72|2|9|8 72|4|9|8
			 96|1|16|6 96|2|16|6 96|3|16|6
			 128|1|16|8 128|2|16|8
			 144|1|16|9 144|2|16|9
			 192|1|16|12
			 256|1|16|16
			 288|1|18|16"
	export BESTCONF=""
	export HPLNB="336"
else
	echo "Unsupported host"
	exit
fi

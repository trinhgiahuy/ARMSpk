#!/bin/bash

export APPDIR="./CoMD"
export BINARY="./bin/CoMD-openmp-mpi"
export INPUT="-iPX -jPY -kPZ"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|96|1|1|1 1|48|1|1|1 1|24|1|1|1 1|12|1|1|1
			 2|24|2|1|1 2|12|2|1|1
			 4|12|2|2|1 4|6|2|2|1
			 12|2|3|2|2
			 24|1|4|3|2
			 32|1|4|4|2
			 48|1|4|4|3
			 96|1|6|4|4"
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

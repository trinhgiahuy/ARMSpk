#!/bin/bash

export APPDIR="./Laghos"
export BINARY="./laghos"
export INPUT="-p 1 -m data/square01_quad.mesh -rs 3 -tf 0.8 -no-vis -pa"
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="6|1 6|2 6|4 8|1 8|2 12|1 12|2 16|1 16|2 24|1 32|1 48|1"
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

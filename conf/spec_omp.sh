#!/bin/bash

export APPDIR="./SPEC_OMP"
#options test->train->ref; 'ref' run for >10min which is too much for gem5
export BINARY="350.md|intel|train
351.bwaves|intel|train
352.nab|intel|train
357.bt331|intel|train
358.botsalgn|intel|train
359.botsspar|intel|train
360.ilbdc|intel|train
362.fma3d|intel|train
363.swim|intel|train
367.imagick|intel|train
370.mgrid331|intel|train
371.applu331|intel|train
372.smithwa|intel|train
376.kdtree|intel|train"
export INPUT=""
export NumRunsTEST=1
export NumRunsBEST=10
export MAXTIME="20m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	export TESTCONF="48"
	export BESTCONF="48"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

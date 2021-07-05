#!/bin/bash

export APPDIR="./SPEC_OMP"
#options test->train->ref; 'ref' run for >10min which is too much for gem5
export BINARY="350.md|train
351.bwaves|train
352.nab|train
357.bt331|train
358.botsalgn|train
359.botsspar|train
360.ilbdc|train
362.fma3d|train
363.swim|train
367.imagick|train
370.mgrid331|train
371.applu331|train
372.smithwa|train
376.kdtree|train"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF=""
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|4 1|8 1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF=""
fi

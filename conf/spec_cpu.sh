#!/bin/bash

export APPDIR="./SPEC_CPU"
#options test->train->ref; 'ref' run for >7min which is too much for gem5
export BINARY="600.perlbench_s|intel|train
602.gcc_s|intel|train
605.mcf_s|intel|train
620.omnetpp_s|intel|train
623.xalancbmk_s|intel|train
625.x264_s|intel|train
631.deepsjeng_s|intel|train
641.leela_s|intel|train
648.exchange2_s|intel|train
657.xz_s|intel|train
603.bwaves_s|intel|train
607.cactuBSSN_s|intel|train
619.lbm_s|intel|train
621.wrf_s|intel|train
627.cam4_s|intel|train
628.pop2_s|intel|train
638.imagick_s|intel|train
644.nab_s|intel|train
649.fotonik3d_s|intel|train
654.roms_s|intel|train"
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

#!/bin/bash

export APPDIR="./SPEC_CPU"
#options test->train->ref; 'ref' run for >7min which is too much for gem5
export BINARY="600.perlbench_s|train
602.gcc_s|train
605.mcf_s|train
620.omnetpp_s|train
623.xalancbmk_s|train
625.x264_s|train
631.deepsjeng_s|train
641.leela_s|train
648.exchange2_s|train
657.xz_s|train
603.bwaves_s|train
607.cactuBSSN_s|train
619.lbm_s|train
621.wrf_s|train
627.cam4_s|train
628.pop2_s|train
638.imagick_s|train
644.nab_s|train
649.fotonik3d_s|train
654.roms_s|train"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="60m"
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
	export TESTCONF="1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF=""
fi

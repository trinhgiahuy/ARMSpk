#!/bin/bash

export APPDIR="./fs2020"
export BINARYS="01.MD.LOOP.makino170209_Ccode.tune03_7loop_oclunroll.170316/main
02.gravity-512.simple.0.N2048_seq_b1_prolog.170316/test
03.kernel_June1_Single_Tune_Full_pr64/kernel_June_1|12
04.kernel_July_Single_c_m1_rvl/kernel_july|12
05.PairList_June_Single_simu2_c_rvl_align/kernel_pairlist_june|12
06.PairList_July_Single_tune_20170901/kernel_pairlist_july|12
08.NICAM.divdamp.r4_collapsed.tune01.170420/postK_nicam_divdamp_ijsplit_tune03
09.NICAM.diffusion.tune_r4.asis_opttune01.170220/postK_nicam_diffusion
10.NICAM.Horizontal_Adv_flux.r4_collapsed.tune03_kpara.30.170418/main
11.NICAM.Horizontal_Adv_limiter.tune_r4.tune01.170407/main
12.NICAM.Vertical_Adv_limiter.tune_r4.asis_opttune01.160902/main
13.NICAM.vi_rhow_solver.tune_r4.asis.160902/main
14.streamlike_pattern1/main
15.streamlike_pattern2/main
16.streamlike_pattern3/main
17.Adventure.region0.tune4pad-acle.armtest.20170821_resize/adventure_kernel_region0_tune4_arm_pad-acle
18.Adventure.region1.tune01.outputcheck.160808.M24/test_dd_dot_product
19.Adventure.region2.tune161011.armtest-M24.170727/test_dd_dot_product|24
20.FFB.callap_kernel2.nodebase_160127.asis.160805/main
21.FFB.spmmv_vec8.170330.funroll_265225.170509/main
22.QCD.ddd_in_s_.qws-0.1.7.tune.170621/main
23.QCD.jinv_ddd_in_s_.qws-0.1.7.tune.170720/main"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="20m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|1 1|4 1|8 1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF="1|12"
	export SCALCONF=""
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF=""
	export BESTCONF=""
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|1 1|4 1|8 1|12 1|16 1|24 1|32 1|36 1|48"
	export BESTCONF="1|12"
elif [ -n "${GEM5HOST}" ]; then
	export GEM5CONF="1|12"
	export NumRunsGEM5=1
fi

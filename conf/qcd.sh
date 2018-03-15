#!/bin/bash

export APPDIR="./QCD/src"
export BINARYS="./ccs_qcd_solver_bench_class1 ./ccs_qcd_solver_bench_class2 ./ccs_qcd_solver_bench_class3"
export BBINARY="./ccs_qcd_solver_bench_class2"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="1m"
export RUNVTUNE="no"

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	# all classes low perf (<10%) when used w/ MPI, fokus: MPI=1 and OMP>1
	export TESTCONF="1|6 1|12 1|24 1|32 1|48 1|96
			 2|12
			 4|6
			 6|4
			 12|2
			 24|1
			 32|1
			 48|1"
	export BESTCONF="1|24"
elif [[ $HOSTNAME = *"lyon"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64 1|128 1|192 1|256
			 4|16
			 16|4
			 32|2
			 64|1
			 96|1
			 128|1
			 192|1
			 256|1"
	export BESTCONF="1|64"
elif [[ $HOSTNAME = *"mill"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

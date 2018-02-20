#!/bin/bash

export APPDIR="./QCD/src"
export BINARYS="./ccs_qcd_solver_bench_class1 ./ccs_qcd_solver_bench_class2 ./ccs_qcd_solver_bench_class3"
export BBINARY="./ccs_qcd_solver_bench_class3"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10

if [[ $HOSTNAME = *"kiev"* ]]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6"
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

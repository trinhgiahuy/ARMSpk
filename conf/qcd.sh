#!/bin/bash

export APPDIR="./QCD/src"
export BINARYS="./ccs_qcd_solver_bench_class1 ./ccs_qcd_solver_bench_class2 ./ccs_qcd_solver_bench_class3"
export INPUT=""

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96 1|48 1|32 1|24 1|12 1|6"
	export BESTCONF=""
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

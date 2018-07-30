#!/bin/bash

export APPDIR="./QCD/src"
export BINARY="./ccs_qcd_solver_bench_class2"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="1m"
export RUNSDE="yes"
export RUNPCM="yes"
export RUNVTUNE="no" #"yes"

if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	# on "normal" Xeon
	# all classes low perf (<10%) when used w/ MPI, fokus: MPI=1 and OMP>1
	export TESTCONF="1|6|1|1|1 1|12|1|1|1 1|24|1|1|1 1|32|1|1|1 1|48|1|1|1 1|96|1|1|1
			 2|6|2|1|1 2|12|2|1|1 2|24|2|1|1
			 4|1|2|2|1 4|2|2|2|1 4|4|2|2|1 4|6|2|2|1 4|12|2|2|1
			 8|1|2|2|2 8|2|2|2|2 8|4|2|2|2 8|6|2|2|2
			 16|1|4|2|2 16|2|4|2|2 16|4|4|2|2
			 24|1|4|3|2 24|2|4|3|2
			 32|1|4|4|2 32|2|4|4|2
			 48|1|4|4|3 48|2|4|4|3
			 64|1|4|4|4"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	# on one of the Phi (knl)
	export TESTCONF="1|64|1|1|1 1|128|1|1|1 1|192|1|1|1 1|256|1|1|1
			 4|16|2|2|1 4|32|2|2|1 4|48|2|2|1 4|64|2|2|1
			 16|4|4|2|2 16|8|4|2|2 16|12|4|2|2 16|16|4|2|2
			 32|2|4|4|2 32|4|4|4|2 32|6|4|4|2 32|8|4|4|2
			 64|1|4|4|4 64|2|4|4|4 64|4|6|4|4 64|6|6|4|4
			 128|1|8|4|4 128|2|8|4|4 128|3|8|4|4
			 256|1|8|8|4 256|2|8|8|4"
	export BESTCONF=""
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	# on one of the Phi (knm)
	export TESTCONF=""
	export BESTCONF=""
else
	echo "Unsupported host"
	exit
fi

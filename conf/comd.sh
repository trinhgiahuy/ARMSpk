#!/bin/bash

export APPDIR="./CoMD"
export BINARY="./bin/CoMD-openmp-mpi"
export INPUT="-iPX -jPY -kPZ" #for 1 MPI rank

if [ "x`lscpu | grep '^Model name.*E5-2650' | wc -l`" = "x1" ]; then
	# on "normal" Xeon
	export TESTCONF="1|96|1|1|1 1|48|1|1|1 1|24|1|1|1 1|12|1|1|1
			 2|24|2|1|1 2|12|2|1|1
			 4|12|2|2|1 4|6|2|2|1
			 12|2|3|2|2
			 24|1|4|3|2
			 32|1|4|4|2
			 48|1|4|4|3
			 96|1|6|4|4"
	export BESTCONF="48|1|4|4|3"
else
	# on one of the Phi
	export TESTCONF=""
	export BESTCONF=""
fi

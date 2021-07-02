#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ SWFFFT =========================================
source conf/swfft.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/swfft/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	# test if Decomposition is valid
	INSIZE="`echo $INPUT | awk '{print $2}'`"
	#XXX: always true for 1==$NumMPI: `dirname $BINARY`/CheckDecomposition $INSIZE $INSIZE $INSIZE $NumMPI > /dev/null 2>&1
	if ! [ "x$?" = "x0" ]; then
		echo "INVALID Decomposition: $BINARY $INPUT" >> $LOG 2>&1
		continue
	fi
	run_gem5_cmd "$@"
done
cd $ROOTDIR

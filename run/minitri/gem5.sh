#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ miniTri ========================================
source conf/minitri.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/minitri/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	if [ "x$NumMPI" != "x1" ]; then
		BINARY=$BINARYMPI
		INPUT=$INPUTMPI
	else
		BINARY=$BINARYOMP
		INPUT="`echo $INPUTOMP | sed -e \"s/OMPNT/$NumOMP/\"`"
	fi
	run_gem5_cmd "$@"
done
cd $ROOTDIR

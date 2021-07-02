#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ XSBench ========================================
source conf/xsbench.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/xsbench/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
DEFINPUT=$INPUT
for BEST in $BESTCONF; do
	NumMPI=1
	INPUT="`echo $DEFINPUT | sed -e \"s/OMPNT/$NumOMP/\"`"
	run_gem5_cmd "$@"
done
cd $ROOTDIR

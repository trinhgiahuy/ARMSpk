#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ miniFE =========================================
source conf/minife.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/minife/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BBINARY; do
		NumMPI=1
		run_gem5_cmd "$@"
		cat ./miniFE.*.yaml >> $LOG 2>&1
		rm -f ./miniFE.*.yaml
	done
done
cd $ROOTDIR

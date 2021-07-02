#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ Laghos =========================================
source conf/laghos.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/laghos/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	run_gem5_cmd "$@"
done
cd $ROOTDIR

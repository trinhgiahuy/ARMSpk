#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ CoMD ===========================================
source conf/comd.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/gem5run/comd/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	X=1; Y=1; Z=1;  #XXX: need to change this if we set NumMPI=1
	#X="`echo $BEST | cut -d '|' -f3`"
	#Y="`echo $BEST | cut -d '|' -f4`"
	#Z="`echo $BEST | cut -d '|' -f5`"
	INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
	run_gem5_cmd "$@"
done
rm *.yaml
cd $ROOTDIR

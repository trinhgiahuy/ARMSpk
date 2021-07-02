#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ FFVC ===========================================
source conf/ffvc.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/gem5run/ffvc/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	X=1; Y=1; Z=1;  #XXX: need to change this if we set NumMPI=1
	#X="`echo $BEST | cut -d '|' -f3`"
	#Y="`echo $BEST | cut -d '|' -f4`"
	#Z="`echo $BEST | cut -d '|' -f5`"
	INPUT=$DEFINPUT
	INPUT="`echo $INPUT | sed -e \"s/DX/$X/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DY/$Y/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DZ/$Z/\"`"
	run_gem5_cmd "$@"
done
cd $ROOTDIR

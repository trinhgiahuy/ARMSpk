#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

export HOSTNAME="kiev0" # we don't care where gem runs
GEM5="$ROOTDIR/dep/gem5_riken/build/ARM/gem5.opt"
GEM5SE="$ROOTDIR/dep/gem5_riken/configs/example/se.py"

# ============================ FFB ============================================
source conf/ffb.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/gem5run/ffb/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	X=1; Y=1; Z=1;  #XXX: need to change this if we set NumMPI=1
	#X="`echo $BEST | cut -d '|' -f3`"
	#Y="`echo $BEST | cut -d '|' -f4`"
	#Z="`echo $BEST | cut -d '|' -f5`"
	INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
	# try finding closest cube size per proc
	FLOAT=`echo "e((1/3)*l($MAXDCZ / $NumMPI))" | bc -l`
	DCZ=`echo "($FLOAT+0.5)/1" | bc`
	INPUT="`echo $INPUT | sed -e \"s/DCZ/$DCZ/\"`"
	mkdir ./tmp${1}; sleep 1; cd ./tmp${1}
	BINARY=../$BINARY
	run_gem5_cmd "$@"
	cd ../; rm -rf ./tmp${1}; sleep 1
done
cd $ROOTDIR

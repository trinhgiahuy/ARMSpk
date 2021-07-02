#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ HPL ============================================
source conf/hpl.sh
LOG="$ROOTDIR/log/`hostname -s`/gem5run/hpl/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
if [ ! -f ./HPL.dat.bak ]; then cp ./HPL.dat ./HPL.dat.bak; fi
for BEST in $BESTCONF; do
	NumMPI=1
	MPIP=1; MPIQ=1; #XXX: always 1 for 1==NumMPI
	#MPIP="`echo $BEST | cut -d '|' -f3`"
	#MPIQ="`echo $BEST | cut -d '|' -f4`"
	sed -e "s/PNS/$HPLNS/" -e "s/PNB/$HPLNB/" -e "s/PPP/$MPIP/" -e "s/PPQ/$MPIQ/" ./HPL.dat.bak > ./HPL.dat
	run_gem5_cmd "$@"
	cat ./HPL.out >> $LOG 2>&1; rm ./HPL.out
done
cd $ROOTDIR

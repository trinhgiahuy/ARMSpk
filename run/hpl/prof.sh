#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please intel-sde-external-8.16.0-2018-01-30-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
if [[ $HOSTNAME = *"kiev"* ]]; then
	SDE="$SDE -bdw -- "
elif [[ $HOSTNAME = *"lyon"* ]]; then
	SDE="$SDE -knl -- "
elif [[ $HOSTNAME = *"mill"* ]]; then
	SDE="$SDE -knm -- "
else
	echo "Unsupported host"
	exit
fi

# ============================ HPL ============================================
source conf/hpl.sh
LOG="$ROOTDIR/log/`hostname -s`/profrun/hpl.log"
mkdir -p `dirname $LOG`
cd $APPDIR
if [ ! -f ./HPL.dat.bak ]; then cp ./HPL.dat ./HPL.dat.bak; fi
for BEST in $BESTCONF; do
	mkdir -p ./oSDE
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	MPIP="`echo $BEST | cut -d '|' -f3`"
	MPIQ="`echo $BEST | cut -d '|' -f4`"
	sed -e "s/PNS/$HPLNS/" -e "s/PNB/$HPLNB/" -e "s/PPP/$MPIP/" -e "s/PPQ/$MPIQ/" ./HPL.dat.bak > ./HPL.dat
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
        mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
	cat ./HPL.out >> $LOG 2>&1; rm ./HPL.out
	for P in `seq 0 $((NumMPI - 1))`; do
                echo "SDE output of MPI process $P" >> $LOG 2>&1
                cat ./oSDE/${P}.txt >> $LOG 2>&1
        done
        echo "=== SDE summary ===" >> $LOG 2>&1
        $ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed -e 's#profrun#bestrun#g'` >> $LOG 2>&1
        rm -rf ./oSDE
done
cd $ROOTDIR
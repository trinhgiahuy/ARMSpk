#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source `cat $ROOTDIR/conf/intel.cfg` intel64 > /dev/null 2>&1
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

# ============================ NICam ==========================================
source conf/nicam.sh
LOG="$ROOTDIR/log/`hostname -s`/profrun/nicam.log"
mkdir -p `dirname $LOG`
cd $APPDIR

# scale down #steps from 11 days to 1 day, and create input data set
sed -i -e 's/^LSMAX  = 0/LSMAX  = 72/'  ../../test.conf
make jobshell > /dev/null 2>&1
if [ -d "$INPUT" ] && [ "x" != "x$INPUT" ]; then cd $INPUT; else exit; fi
export FORT_FMT_RECL=400
ln -s ../../../../bin/nhm_driver .
ln -s ../../../../data/mnginfo/rl00-prc10.info .
ln -s ../../../../data/grid/vgrid/vgrid40_24000-600m.dat .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000000 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000001 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000002 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000003 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000004 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000005 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000006 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000007 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000008 .
ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000009 .

for BEST in $BESTCONF; do
	mkdir -p ./oSDE
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY\"" >> $LOG 2>&1
	mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY" >> $LOG 2>&1
	cat ./msg.pe00000 >> $LOG 2>&1
	for P in `seq 0 $((NumMPI - 1))`; do
		echo "SDE output of MPI process $P" >> $LOG 2>&1
		cat ./oSDE/${P}.txt >> $LOG 2>&1
	done
	echo "=== SDE summary ===" >> $LOG 2>&1
	$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed 's#profrun#bestrun#g'` >> $LOG 2>&1
	rm -rf ./oSDE
	if [ "x$RUNVTUNE" = "xyes" ]; then
		echo "=== vtune hpc-performance ===" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY
		amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
		rm -rf ./oVTP.`hostname`
		echo "=== vtune memory-access ===" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY
		amplxe-cl -report summary -q -result-dir ./oVTM.`hostname` >> $LOG 2>&1
		rm -rf ./oVTM.`hostname`
	fi
done

# clean up
cd $ROOTDIR
cd $APPDIR
if [ -d "$INPUT" ] && [ "x" != "x$INPUT" ]; then rm -rf $INPUT; fi

cd $ROOTDIR

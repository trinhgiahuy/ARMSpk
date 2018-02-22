#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
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
LOG="$ROOTDIR/log/bestrun/nicam.log"
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
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		cat ./msg.pe00000 >> $LOG 2>&1
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done

# clean up
cd $ROOTDIR
cd $APPDIR
if [ -d "$INPUT" ] && [ "x" != "x$INPUT" ]; then rm -rf $INPUT; fi

echo "Best NICam run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR

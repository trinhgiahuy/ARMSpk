#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
export PATH=$ROOTDIR/dep/likwid/bin:$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then UNCORE="--umin 2.7 --umax 2.7"; fi

# ============================ NICam ==========================================
source conf/nicam.sh
LOG="$ROOTDIR/log/`hostname -s`/freqrun/nicam.log"
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

for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
	for BEST in $BESTCONF; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			cat ./msg.pe00000 >> $LOG 2>&1
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done

# clean up
cd $ROOTDIR
cd $APPDIR
if [ -d "$INPUT" ] && [ "x" != "x$INPUT" ]; then rm -rf $INPUT; fi

cd $ROOTDIR

#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
#source $ROOTDIR/conf/intel.cfg
#source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
#source $ROOTDIR/dep/spack/share/spack/setup-env.sh
#spack load gcc@8.4.0
### all static, no need for compilers and envs

SPECCMD="runcpu --config=nedo.cfg --nobuild --action=run --noreportable --use_submit_for_speed"

SDEPATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin
if [ ! -x "`PATH=$SDEPATH:$PATH which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;

# ============================ SPEC CPU =======================================
source conf/spec_cpu.sh
LOGDIR="$ROOTDIR/log/`hostname -s`/profrun/spec_cpu"
mkdir -p $LOGDIR
cd $APPDIR
for BENCH in $BINARY; do
	BM="`echo $BENCH | cut -d '|' -f1`"
	COMP="`echo $BENCH | cut -d '|' -f2`"
	SIZE="`echo $BENCH | cut -d '|' -f3`"
	LOG="$LOGDIR/${BM}.log"
	NumMPI=1
	for NumOMP in $BESTCONF; do
		if [ "x$RUNSDE" = "xyes" ]; then
			echo "=== sde run for $BENCH ===" >> $LOG 2>&1
			mkdir -p $LOGDIR/${BM}_${NumMPI}_${NumOMP}_sde
			START="`date +%s.%N`"
			bash -c "export PATH=$SDEPATH:$PATH; source ./shrc; $SPECCMD --iterations=1 --size=$SIZE --threads=$NumOMP --define COMP=sde --define RESDIR=$LOGDIR/${BM}_${NumMPI}_${NumOMP}_sde $BM" >> $LOG 2>&1
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			# in case the RESDIR stuff fails
			mkdir -p $LOGDIR/${BM}_${NumMPI}_${NumOMP}_sde/extra
			find . -name 'dcfg-out.*' -exec mv {} $LOGDIR/${BM}_${NumMPI}_${NumOMP}_sde/extra \;
			REPORT="`find $ROOTDIR/$APPDIR/result -type f -name '*.log' | sort -g | tail -1`"
			if ! grep '^Success:' $LOG > /dev/null 2>&1 ; then echo "=== error report ===" >> $LOG 2>&1; cat $REPORT >> $LOG 2>&1; fi
		fi
	done
done
cd $ROOTDIR

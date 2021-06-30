#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

# ============================ FS2020 ====================================
source conf/fs2020.sh
LOGDIR="$ROOTDIR/log/`hostname -s`/bestrun/fs2020"
mkdir -p $LOGDIR
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	NumOMP=$BEST
	moreMPI="-x OMP_NUM_PARALELL=$NumOMP -x OMP_PROC_BIND=close -x FLIB_FASTOMP=FALSE -x FLIB_CNTL_BARRIER_ERR=FALSE"
	for BMconf in $BINARYS; do
		BINARY="`echo ${BMconf} | cut -d '|' -f1`"
		MAXOMP="`echo ${BMconf} | cut -s -d '|' -f2`"
		if [ ! -z $MAXOMP ] && [ $NumOMP -gt $MAXOMP ]; then continue; fi
		if [[ "$BINARY" = "23."* ]] && [ $NumOMP -lt 6 ]; then continue; fi	# needs at least 6 threads to avoid floating-point exception
		LOG="${LOGDIR}/`echo ${BINARY} | cut -d'/' -f1`.${NumOMP}.log"
		echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG $moreMPI) $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG $moreMPI) $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
		BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sed -e 's/D/e/g' | sort -g | head -1`"
		echo "Best $BINARY NumOMP=$NumOMP run: $BEST"
	done
done
cd $ROOTDIR

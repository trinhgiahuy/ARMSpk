#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

# ============================ BabelStream ====================================
source conf/babelstream.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/bestrun/babelstream"
mkdir -p `dirname $DEFLOG`
cd $APPDIR
DEFINPUT=$INPUT
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		NumOMP=$BEST
		S="`echo $BINARY | cut -d '_' -f2`"
		BINARY="`echo $BINARY | cut -d '_' -f1`"
		LOG="${DEFLOG}${S}gb.log"
		S=$((S*1024*1024*1024/8))
		INPUT="`echo $DEFINPUT | sed -e \"s/SIZE/$S/\"`"
		echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done
echo "Best BabelStream run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR

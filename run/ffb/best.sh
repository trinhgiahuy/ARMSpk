#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ FFB ============================================
source conf/ffb.sh
NumRUNS=10
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/bestrun/ffb.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	X="`echo $BEST | cut -d '|' -f3`"
	Y="`echo $BEST | cut -d '|' -f4`"
	Z="`echo $BEST | cut -d '|' -f5`"
	INPUT=$DEFINPUT
	INPUT="`echo $INPUT | sed -e \"s/DX/$X/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DY/$Y/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DZ/$Z/\"`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mkdir ./tmp; sleep 2; cd ./tmp
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1 
		cd ../; rm -rf ./tmp; sleep 2
	done
done
echo "Best FFB run:"
BEST="`grep '^MAIN_LOOP' $LOG | awk -F 'max)' '{print $2}' | awk -F '[' '{print $1}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
